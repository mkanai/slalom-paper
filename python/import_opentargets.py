import argparse
import hail as hl
from hail.utils.java import Env
from fm_insights.utils import register_log

bucket = "finucane-opentargets-genetics"

# gsutil -u encode-uk-biobank-restrict -m rsync -rd gs://open-targets-genetics-releases/22.02.01/variant-index gs://finucane-opentargets-genetics/releases/22.02.01/variant-index

def main(args):
    if args.import_spark:
        spark = Env.spark_session()
        datasets = ["lut/study-index", "lut/variant-index", "v2d", "v2d_credset", "variant-index"]

        for dataset in datasets:
            _reader = spark.read.parquet if dataset == "variant-index" else spark.read.json
            df = _reader(f"gs://{bucket}/releases/{args.release}/{dataset}/")
            ht = hl.Table.from_spark(df)
            ht = ht.repartition(1000)
            ht = ht.checkpoint(f"gs://{bucket}/releases/{args.release}/{dataset}.ht", overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=str, default="21.10.01")
    parser.add_argument("--import-spark", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
