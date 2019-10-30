
def index(args):
    """create an index for SEQC.

    :param args: parsed arguments. This function is only called if subprocess_name is
      'index'
    """

    # functions to be pickled and run remotely must import all their own modules
    import sys
    import logging
    from seqc import ec2, log
    from seqc.sequence.index import Index
    from seqc import version

    logging.basicConfig(
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(args.log_name),
            logging.StreamHandler(sys.stdout)
        ]
    )

    log.info("SEQC v{}".format(version.__version__))
    log.args(args)

    idx = Index(args.organism, args.ids)
    idx.create_index(
        s3_location=args.upload_prefix,
        ensemble_release=args.ensemble_release,
        read_length=args.read_length
    )

    log.info("DONE.")
