
def index(args):
    """create an index for SEQC.

    :param args: parsed arguments. This function is only called if subprocess_name is
      'index'
    """

    # functions to be pickled and run remotely must import all their own modules
    from seqc import ec2, log
    from seqc.sequence.index import Index

    log.setup_logger(args.log_name)
    log.args(args)

    if args.remote:
        with ec2.instance_clean_up(args.email, args.upload_prefix, log_name=args.log_name):
            idx = Index(args.organism, args.ids)
            idx.create_index(
                s3_location=args.upload_prefix,
                ensemble_release=args.ensemble_release,
                read_length=args.read_length
            )
    else:
        idx = Index(args.organism, args.ids)
        idx.create_index(
            s3_location=args.upload_prefix,
            ensemble_release=args.ensemble_release,
            read_length=args.read_length
        )
