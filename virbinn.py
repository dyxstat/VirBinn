#########The structure of the main script is modified from bin3C########
from Script.raw_contact import ContactMatrix
from Script.impute import impute
from Script.cluster import bin
from Script.utils import make_dir
from Script.exceptions import ApplicationException
import argparse
import warnings
import logging
import shutil
import sys
import os


##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")

__version__ = '1.0.0, released at 05/2025'

if __name__ == '__main__':
    
    def mk_version():
        return 'VirBinn v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg

    runtime_defaults = {
        'min_len': 1000,
        'min_signal': 1,
        'min_mapq': 30,
        'min_match': 30,
        'discard_viral': 50,
        'discard_host': 80
    }

    script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/VirBinn.log]')


    parser = argparse.ArgumentParser(description='VirBinn: binning viral contigs through imputation')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')

    cmd_raw = subparsers.add_parser('raw', parents=[global_parser],
                                      description='Raw contacts.')
       
    cmd_impute = subparsers.add_parser('impute', parents=[global_parser],
                                      description='Do the imputation.')
                                   
    cmd_cl = subparsers.add_parser('bin', parents=[global_parser],
                                      description='Do the binning.')

    cmd_test = subparsers.add_parser('test', parents=[global_parser],
                                        description='pipeline testing.')

    '''
    Generating raw contacts subparser input
    '''
    cmd_raw.add_argument('--min-len', type=int,
                           help='Minimum acceptable contig length [1000]')
    cmd_raw.add_argument('--min-signal', type=int,
                           help='Minimum acceptable Hi-C signal [1]')
    cmd_raw.add_argument('--min-mapq', type=int,
                           help='Minimum acceptable mapping quality [30]')
    cmd_raw.add_argument('--min-match', type=int,
                           help='Accepted alignments must being N matches [30]')
    cmd_raw.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
                           help='Case-sensitive enzyme name. Use multiple times for multiple enzymes')
    cmd_raw.add_argument('FASTA', help='Reference fasta sequence')
    cmd_raw.add_argument('BAM', help='Input bam file in query order')
    cmd_raw.add_argument('OUTDIR', help='Output directory')
    
         
    '''
    Do the imputation
    '''
    cmd_impute.add_argument('--discard-viral', type=float, default=None, help='Retention ratio for edges')
    cmd_impute.add_argument('--discard-host', type=float, default=None, help='Retention ratio for edges')
    cmd_impute.add_argument('VIRAL', help='viral contig name.txt')
    cmd_impute.add_argument('OUTDIR', help='Output directory')
    
       
    '''
    Clutering subsparser input
    '''
    cmd_cl.add_argument('--precompute', default=False, action='store_true', help='use precomputed matrix to compute score')
    cmd_cl.add_argument('--output-prefix', type=str, default='viral2bin',
                               help='output prefix')
    cmd_cl.add_argument('FASTA', help='Reference fasta sequence')
    cmd_cl.add_argument('OUTDIR', help='Output directory of sub bins')
    
    
    '''
    Testing of the VirBinn software
    '''
    cmd_test.add_argument('--OUTDIR', type=str, default='Test/out_test', help='Output directory of testing results')

    args = parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.cover)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)
    
    # Create temp folder
    temp_folder = os.path.join(args.OUTDIR , 'tmp')
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)
    else:
        shutil.rmtree(temp_folder)           
        os.mkdir(temp_folder)
           
    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join(args.OUTDIR, 'VirBinn.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:
        if args.command == 'raw':
            if args.enzyme is not None:
                logger.info('Begin constructing raw contact matrix...')
                cm = ContactMatrix(args.BAM,
                                args.enzyme,
                                args.FASTA,
                                args.OUTDIR,
                                min_mapq=ifelse(args.min_mapq, runtime_defaults['min_mapq']),
                                min_len=ifelse(args.min_len, runtime_defaults['min_len']),
                                min_match=ifelse(args.min_match, runtime_defaults['min_match']),
                                min_signal=ifelse(args.min_signal, runtime_defaults['min_signal']))

                logger.info('Raw contact matrix construction finished')
                
        elif args.command == 'impute':
            logger.info('Begin imputation stage...')
            impute(args.VIRAL,
                   args.OUTDIR,
                   discard_viral=ifelse(args.discard_viral, runtime_defaults['discard_viral']),
                   discard_host=ifelse(args.discard_host, runtime_defaults['discard_host']))
                   
            logger.info('Imputation stage finished')
            
            
        elif args.command == 'bin':
            logger.info('Begin binning stage...')
            bin(args.FASTA, 
                args.OUTDIR, 
                precompute = args.precompute, 
                output_prefix = 'viral2bin')
            
            bin(fasta=args.FASTA,
                       cluster_file=args.cluster_input,
                       output_prefix=args.output_prefix,
                       outdir=args.OUTDIR)
            logger.info('Binning stage finished')



    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)
