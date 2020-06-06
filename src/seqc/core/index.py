
def index(args):
    """create an index for SEQC.

    :param args: parsed arguments. This function is only called if subprocess_name is
      'index'
    """

    # functions to be pickled and run remotely must import all their own modules
    import sys
    import logging
    from seqc import ec2, log, io
    from seqc.sequence.index import Index
    from seqc.alignment import star
    from seqc import version

    logging.basicConfig(
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(args.log_name),
            logging.StreamHandler(sys.stdout)
        ]
    )

    log.info("SEQC=v{}".format(version.__version__))
    log.info("STAR=v{}".format(star.get_version()))
    log.args(args)

    with ec2.instance_clean_up(
        email=args.email, upload=args.upload_prefix, log_name=args.log_name,
        debug=args.debug, terminate=args.terminate, running_remote=args.remote
    ):

        idx = Index(args.organism, args.ids, args.folder)
        idx.create_index(
            s3_location=args.upload_prefix,
            ensemble_release=args.ensemble_release,
            read_length=args.read_length,
            valid_biotypes=args.valid_biotypes
        )

        # upload the log file (seqc_log.txt, nohup.log, Log.out)
        if args.upload_prefix:
            bucket, key = io.S3.split_link(args.upload_prefix)
            for item in [args.log_name, "./nohup.log", "./Log.out"]:
                try:
                    ec2.Retry(retries=5)(io.S3.upload_file)(
                        item, bucket, key
                    )
                    log.info("Successfully uploaded {} to the specified S3 location {}".format(item, args.upload_prefix))
                except FileNotFoundError:
                    log.notify("Item {} was not found! Continuing with upload...".format(item))

    log.info("DONE.")
