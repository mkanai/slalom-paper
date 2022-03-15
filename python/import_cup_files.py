import hail as hl

for reference_genome in ["GRCh37", "GRCh38"]:
    for file in hl.hadoop_ls(f"gs://meta-finemapping-simulation/cup_files/*{reference_genome}.*.bed"):
        bed_path = file["path"]
        ht_path = bed_path.replace(".bed", ".ht")
        bed = hl.import_bed(bed_path, reference_genome=reference_genome, min_partitions=10)
        bed.write(ht_path)
