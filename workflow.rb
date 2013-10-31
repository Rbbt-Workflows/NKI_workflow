require 'rbbt/workflow'
require 'rbbt/entity/study'
require 'rbbt/entity/study/genotypes'

module NKIWorkflow
  extend Workflow

  def self.exome_bed_file_for_genes(genes)
    bed = TSV.setup({}, :key_field => "Ensembl Exon ID", :fields => %w(chr start end gene), :type => :list)
    organism = genes.organism
    exon_info = Organism.exons(organism).tsv :fields => ["Chromosome Name", "Exon Chr Start", "Exon Chr End"], :persist => true
    exons_for = Organism.exons(organism).tsv :key_field => "Ensembl Gene ID", :fields => ["Ensembl Exon ID"], :persist => true, :merge => true, :type => :flat
    genes.each do |gene|
      exons = exons_for[gene.ensembl]
      exons.each do |exon|
        chr, start, eend = exon_info[exon]
        bed[exon] = [chr, start, eend, gene]
      end
    end
    bed
  end

  input :matrix, :tsv, "Gene X Sample binary matrix"
  input :max, :integer, "Max number of genes to consider", 20
  task :exclusivity => :tsv do |matrix, max|
    script =<<-EOF 

    source('#{Rbbt.share.find(:lib).R["MM.assess.mutation.patterns.R"]}'); 
    m = as.matrix(t(data)); 
    m[is.na(m)] = FALSE; 
    data = assess.me.event.matrix(m, #{max});

    data
    EOF

    matrix.R(script)
  end

  input :study, :string, "Study", nil
  input :threshold, :float, "P-value threshold", 0.1
  task :significantly_mutated => :tsv do |study,threshold|
    study = Study.setup(study.dup)

    matrix = study.gene_sample_matrix
    matrix.add_field "Associated Gene Name" do |gene, values|
      gene.name
    end

    matrix = matrix.reorder "Associated Gene Name"
    matrix = matrix.slice(matrix.fields - ["Ensembl Gene ID"])

    bedfile = study.dir.data["kinome-bed.tsv"]
    if not bedfile.exists?
      bedfile = self.file("bedfile")
      bed_tsv = NKIWorkflow.exome_bed_file_for_genes(matrix.keys.name)
      Open.write(bedfile) do |file|
        file.puts bed_tsv.fields * "\t"
        bed_tsv.each do |key, values|
          file.puts [key].concat(values) * "\t"
        end
      end
    end

    script = "source('#{Rbbt.share.find(:lib).R["MM.identify.mutated.genes.binom.R"]}'); m = as.matrix(t(data)); m[is.na(m)] = FALSE; bed = get.sizes.from.BED('#{bedfile}'); data = run.identification.signif.genes(m, bed);"
    tsv = matrix.R(script).select("p.value"){|v| v.to_f < threshold}
    tsv.fields = tsv.fields.collect{|f| f == "gene" ? "Associated Gene Name" : f }
    tsv = tsv.reorder("Associated Gene Name", tsv.fields - ["Associated Gene Name"])
    tsv.namespace = study.organism
    tsv
  end

  export_asynchronous :significantly_mutated
end
