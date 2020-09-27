import argparse
import logging
import os
import json
import simulate_DM

parser = argparse.ArgumentParser(description='DM simulation script.')
parser.add_argument('-o', '--output_dir', help='Output directory', required=True)
parser.add_argument('-p', '--json_params', help='Params in JSON format', required=True)
parser.add_argument('-j', '--job_id', help='Job ID', required=False, default='')
parser.add_argument('-r', '--random_seed', help='Random seed', required=False, default=0)

def generate_logger(job_dir, job_id):
    # Initialise logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create file handler which logs even debug messages
    log_file_path = os.path.join(job_dir, 'log_{0}.txt'.format(job_id))
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(fmt='%(asctime)s [%(levelname)s]: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info('Created log at path {path}.'.format(path=log_file_path))

    return logger

def run():
    args = parser.parse_args()

    # Define job directory
    job_dir = args.output_dir

    # Create job dir if it doesn't exist
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)


    # Initialise logger
    logger = generate_logger(job_dir, job_id = str(args.job_id))

    # Start logging
    logger.info('\n------\nSetting up job {job_id}\n------\n'.format(job_id=args.job_id))
    logger.info('User-defined arguments:')
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # Run simulations
    sim_job = simulate_DM.SimulationJob(path_json_params=args.json_params, outpath=job_dir, logger=logger, job_id=args.job_id, random_seed=args.random_seed)
    sim_job.run()

if __name__ == '__main__':
    run()
