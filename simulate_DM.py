import uuid, json, os
import random
from random import randint, shuffle, sample
import math
from numpy.random import binomial
from numpy import log2
import pandas as pd
import matplotlib.pyplot as plt
import logging
import time
from pprint import pprint, pformat

# Define params
num_cells = 2
num_cycles = 200
num_simulations = 100

param_p_tandup = 0.01

min_cassette = 1
max_cassette = 100

max_randnum = 1000


class Cell(object):

    def __init__(self, DM_lst=[1], parent=None, p_tandup=0.0, strain=None, generation_cycle=None):
        self.id = uuid.uuid4().hex
        self.DM = DM_lst
        self.parent = parent
        self.p_tandup = p_tandup
        self.strain = strain
        self.generation_cycle = generation_cycle

    def divide(self):
        DM_lst = self.DM
        shuffle(DM_lst)

        if len(DM_lst) > 1:
            DM_lst = self.tandem_duplication(DM_lst)

        # Duplicate DMs
        DM_lst = DM_lst * 2

        # Split DMs for daughter cells
        split_idx = binomial(len(DM_lst), 0.5)

        # Construct new cells
        daughter1 = Cell(DM_lst=DM_lst[:split_idx], parent=self.id, p_tandup=self.p_tandup, strain=self.strain)
        daughter2 = Cell(DM_lst=DM_lst[split_idx:], parent=self.id, p_tandup=self.p_tandup, strain=self.strain)

        return daughter1, daughter2

    def tandem_duplication(self, DM_lst):
        if self.p_tandup == 0.0:
            return DM_lst

        # Pick a DM at random
        this_rand = randint(0, int(1/self.p_tandup))
        if this_rand < len(DM_lst):
            # Pick the DM and duplicate a random part of it
            this_DM = DM_lst[this_rand]
            this_DM = this_DM + randint(0, this_DM)
            new_DM_lst = [DM_lst[i] for i in range(0, len(DM_lst)) if i != this_rand]
            new_DM_lst.append(this_DM)
        else:
            new_DM_lst = DM_lst

        return new_DM_lst

    def to_table(self):
        return pd.DataFrame({'id': self.id, 'parent': self.parent, 'DM': self.DM, 'generation_cycle': self.generation_cycle})

    def to_summary(self):
        if len(self.DM) == 0:
            max_cass_count = 0
            total_cass = 0
            num_dm = 0
        else:
            max_cass_count = max(self.DM)
            total_cass = sum(self.DM)
            num_dm = len(self.DM)
        return {'id': self.id, 'parent': self.parent, 'generation_cycle': self.generation_cycle,
                'max_cass_count': max_cass_count, 'total_cass': total_cass, 'num_dm': num_dm}

    def __int__(self):
        return sum(self.DM)

    def __float__(self):
        return float(self.__int__())

    def __add__(self, other):
        return float(self) + float(other)

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __repr__(self):
        return "ID: {0}, parent: {1}, DM number: {2}, DM list: {3}, td probability: {4}".format(self.id, self.parent, sum(self.DM), str(self.DM), str(self.p_tandup))


class Simulation(object):

    def __init__(self, num_cells, param_p_tandup, DM_init=[1], resistance_fx='linear', strains=None, max_num_cells=1000,
                 max_num_cycles=1000, logger=None):
        self.num_cells = num_cells
        self.param_p_tandup = param_p_tandup
        self.DM_init=DM_init
        self.resistance_fx = resistance_fx
        self.max_num_cells = max_num_cells
        self.max_num_cycles = max_num_cycles
        self.logger = logger

        if strains:
            if len(strains) != num_cells:
                raise Exception('Number of strains does not correspond to number of cells.')
            else:
                self.strains = strains
        cells = []

        for i in range(0, num_cells):
            if strains:
                this_strain = strains[i]
            else:
                this_strain = None
            #print DM_init
            cells.append(Cell(p_tandup=param_p_tandup, DM_lst=DM_init, strain=this_strain, generation_cycle=0))
        #print cells
        self.cells = cells

    def run(self):
        # Loop through cycles
        cycle = 0
        num_cells = self.num_cells
        pre_cells = self.cells
        max_num_cells = self.max_num_cells
        max_num_cycles = self.max_num_cycles
        #print self.cells

        run_data = pre_cells
        cells = pre_cells
        cycle_summary = []
        start_time = time.time()
        while cycle < max_num_cycles:
            cycle+=1

            # Purge cells below cassette threshold
            survivors = [cell for cell in cells if float(cell) >= min_cassette]
            cells = survivors
            if self.logger:
                self.logger.debug('Cycle {0}, number of cells: {1}'.format(str(cycle), str(len(cells))))

            # Log every 100th cycle
            if self.logger:
                if cycle%10 == 0:
                    end_time = time.time()
                    self.logger.info('Cycle {0}, number of cells: {1}. Time elapsed: {el:2.3f}s'.format(str(cycle), str(len(cells)), el=end_time-start_time))
                    start_time = end_time

            # Subsample population if larger than num_cells
            pre_cell_num = len(cells) #+ len(run_data)

            if pre_cell_num > max_num_cells:
                cell_excess = pre_cell_num - max_num_cells
                # Generate random subset

                sample_idx = sample(range(0, pre_cell_num), max_num_cells)
                #cells = [cells[i - len(run_data)] for i in sample_idx if i > (len(run_data)-1)]
                cells = [cells[i] for i in sample_idx]
                #run_data = [run_data[i] for i in sample_idx if i < len(run_data)]
                if self.logger:
                    self.logger.debug('Subsetting cell population from {0} to {1} (cycle #{2}).'.format(str(pre_cell_num), str(len(cells)), str(cycle)))
            else:
                cell_excess = 0

            # Initialise random numbers for current cycle
            rand_numbers = [randint(0, max_cassette) for i in range(0, len(cells))]

            # Loop through cells
            post_div_cells = []
            for i in range(0, len(cells)):
                this_cell = cells[i]

                # Decide whether to divide
                division_val = self.get_resistance(this_cell)
                if division_val > rand_numbers[i]:
                    daughter1, daughter2 = this_cell.divide()
                    daughter1.generation_cycle = cycle
                    daughter2.generation_cycle = cycle
                    post_div_cells.append(daughter1)
                    post_div_cells.append(daughter2)
                else:
                    post_div_cells.append(this_cell)

            cells = post_div_cells
            if len(post_div_cells) > 0:
                post_div_tbl = pd.DataFrame([cell.to_summary() for cell in post_div_cells])
                post_div_tbl = post_div_tbl.loc[post_div_tbl['max_cass_count'] > 0]
                this_summary = self.get_generation_summary(post_div_tbl)
                this_summary['cell_excess'] = cell_excess
                cycle_summary.append(this_summary)
            #run_data += post_div_cells

        self.cells = cells #run_data
        self.cycle_summary = cycle_summary

    def get_resistance(self, x):
        if self.resistance_fx == 'linear':
            return float(x)
        if self.resistance_fx == 'asympt':
            return float(x)/(float(x) + 10)
        if self.resistance_fx == 'log':
            return log2(float(x))
        if self.resistance_fx == 'exp':
            return math.exp(float(x))

    def summarise(self):
        if self.logger:
            self.logger.info('Concatenating {n} cells.'.format(n = str(len(self.cells))))
        cell_tbl = pd.concat([cell.to_table() for cell in self.cells])
        if self.logger:
            self.logger.info('Calculating cassette sums.')
        dm_sums = cell_tbl.groupby('id', as_index=False).sum().rename(columns={'DM': 'dm_sums'})
        if self.logger:
            self.logger.info('Calculating cassette averages.')
        dm_means = cell_tbl.groupby('id', as_index=False).mean().rename(columns={'DM': 'dm_means'})
        if self.logger:
            self.logger.info('Calculating DM numbers.')
        dm_num = cell_tbl.groupby('id', as_index=False).count().rename(columns={'DM': 'dm_num'})
        summary_tbl = dm_sums #pd.merge(dm_sums, dm_means, how='left', on='id')
        summary_tbl['p_tandup'] = self.param_p_tandup
        summary_tbl['DM_init'] = str(self.DM_init)
        summary_tbl['resistance_fx'] = self.resistance_fx

        # Create dm_tbl where each DM gets its own row
        dm_summary = self.get_dm_summary(self.cells)
        dm_summary = self.summarise_cycles(dm_summary)

        # Create cycle summary
        cycle_summary = self.summarise_cycles(pd.concat(self.cycle_summary))

        return cell_tbl, summary_tbl, dm_summary, cycle_summary

    def get_dm_summary(self, cells):
        # Generate table from list of cells
        dm_tbl = pd.concat([cell.to_table() for cell in cells])
        dm_tbl.columns = ['cassette_count', 'generation_cycle', 'id', 'parent']
        dm_tbl.reset_index(inplace=True, drop=True)

        # Generate dm_summary table
        dm_summary = dm_tbl

        # Group by cell ID to assign maximum cassette count and total cassette count per cell
        dm_summary = dm_summary.groupby(by=['id', 'generation_cycle']).agg({'cassette_count': ['max', 'sum', 'count']})
        dm_summary = pd.DataFrame(dm_summary.to_records())
        dm_summary.columns = ['id', 'generation_cycle', 'max_cass_count', 'total_cass', 'num_dm']

        # Summarise each generation cycle, ie. get the number of cells and cassettes in each cycle
        #dm_summary = dm_summary.groupby(by=['generation_cycle', 'max_cass_count']).agg(
        #    {'id': 'count', 'total_cass': 'sum'}).reset_index()
        #dm_summary.columns = ['generation_cycle', 'max_cass_count', 'num_cells', 'total_cass']
        dm_summary = self.get_generation_summary(dm_summary)

        return dm_summary

    def get_generation_summary(self, cells):
        # Summarise each generation cycle, ie. get the number of cells and cassettes in each cycle
        generation_summary = cells.groupby(by=['generation_cycle', 'max_cass_count']).agg(
            {'id': 'count', 'total_cass': 'sum', 'num_dm': 'sum'}).reset_index()
        generation_summary.columns = ['generation_cycle', 'max_cass_count', 'num_dm', 'num_cells', 'total_cass']
        return generation_summary

    def summarise_cycles(self, cycle_summary):

        # Sort by generation cycle, assign cumulative cell count
        cycle_groups = cycle_summary.groupby(by='max_cass_count')
        group_lst = []
        for name, group in cycle_groups:
            group = group.sort_values(by='generation_cycle')
            group['count_sum'] = group['num_cells'].cumsum()
            group['cass_sum'] = group['total_cass'].cumsum()
            group['dm_sum'] = group['num_dm'].cumsum()
            group_lst.append(group)
        cycles = pd.concat(group_lst)
        cycles.reset_index(inplace=True, drop=True)

        return cycles

    def plot(self, p1 = 'output/hist.pdf', p2 = 'output/sums_vs_means.pdf'):
        cell_tbl, summary_tbl = self.summarise()

        # Plot histogram
        num_bins = int(max(summary_tbl['dm_sums'])/2)
        if num_bins < 1:
            num_bins = 1
        fig = plt.figure(1)
        plt.hist(summary_tbl['dm_sums'], bins=num_bins)
        fig.savefig(p1)
        plt.close(fig)

        fig = plt.figure(2)
        plt.plot(summary_tbl['dm_sums'], summary_tbl['dm_means'], 'k+')
        plt.xlabel('Total cassettes')
        plt.ylabel('Mean cassettes per DM')
        fig.savefig(p2)
        plt.close(fig)


class SimulationJob(object):

    def __init__(self, path_json_params, outpath, logger, job_id=None, random_seed=0):
        self.path_json_params = path_json_params
        self.outpath = outpath
        self.logger = logger
        if job_id:
            self.job_id = job_id
        else:
            self.job_id = ''
        self.rand_seed = float(random_seed)

        self.logger.info('Initialising simulation job.')

        self.logger.info('Reading parameters from {0}'.format(path_json_params))
        self.set_params_from_json(path_json_params=path_json_params)

        # Set random seed
        random.seed(self.rand_seed)

        # Initialise output files
        self.output_sim_param = os.path.join(self.outpath, 'sim_tbl_multi_'+self.job_id+'.csv')
        self.outpath_dm = os.path.join(self.outpath, 'dm_tbl_multi_'+self.job_id+'.csv')
        self.outpath_cycles = os.path.join(self.outpath, 'cycles_multi_'+self.job_id+'.csv')
        for f in [self.output_sim_param, self.outpath_dm, self.outpath_cycles]:
            open(f, 'w').close()

    def set_params_from_json(self, path_json_params):
        json_params = self.load_json_params(path_json_params)

        global_params = json_params['global']
        self.num_cells = global_params['num_cells']
        self.num_cycles = global_params['num_cycles']
        self.num_simulations = global_params['num_simulations']
        self.param_p_tandup = global_params['param_p_tandup']
        self.min_cassette = global_params['min_cassette']
        self.max_cassette = global_params['max_cassette']
        self.max_randnum = global_params['max_randnum']
        self.max_num_cells = global_params['max_num_cells']

        experiments = json_params['experiments']

        experiments_parsed = []
        for experiment in experiments:
            if experiment['strains'] == 'None':
                experiment['strains'] = None
            experiment['DM_init'] = eval(experiment['DM_init'])

            experiments_parsed.append(experiment)
        self.experiments = experiments_parsed

    def load_json_params(self, path_json_params):
        with open(path_json_params) as handle:
            param_dict = json.loads(handle.read())

        return param_dict

    def run(self):
        num_simulations = self.num_simulations
        experiments = self.experiments

        num_experiments = len(experiments)

        sim_param_lst = []
        dm_tbl_last_allparam = []
        for sim_idx in range(0, len(experiments)):
            self.logger.info('Running parameter setting {0} of {1}'.format(str(sim_idx + 1), str(num_experiments)))
            this_experiment = experiments[sim_idx]
            self.logger.info(this_experiment)
            sim_lst = []
            dm_tbl_lst = []
            for sim in range(0, num_simulations):
                self.logger.info('Running simulation number {0} of {1}'.format(str(sim + 1), str(num_simulations)))
                simulation = Simulation(num_cells=num_cells,
                                        param_p_tandup=this_experiment['p_tandup'],
                                        DM_init=this_experiment['DM_init'],
                                        resistance_fx=this_experiment['resistance_fx'],
                                        strains=this_experiment['strains'],
                                        max_num_cells=self.max_num_cells,
                                        max_num_cycles=self.num_cycles,
                                        logger=self.logger)
                simulation.run()
                cell_tbl, summary_tbl, dm_tbl, cycle_summary = simulation.summarise()
                summary_tbl['dm_sums'].astype('int', inplace=True)
                summary_tbl['sim'] = sim
                summary_tbl['experiment'] = this_experiment['experiment']
                #sim_lst.append(summary_tbl)
                dm_tbl['sim'] = sim
                dm_tbl['experiment'] = this_experiment['experiment']
                dm_tbl['job_id'] = self.job_id
                #dm_tbl_lst.append(dm_tbl)
                cycle_summary['sim'] = sim
                cycle_summary['experiment'] = this_experiment['experiment']
                cycle_summary['job_id'] = self.job_id

                # Append simulation results to output files
                if sim == 0 and sim_idx == 0:
                    param_header = True
                else:
                    param_header = False
                with open(self.output_sim_param, 'a') as f:
                    summary_tbl.to_csv(f, header=param_header)
                with open(self.outpath_dm, 'a') as f:
                    dm_tbl.to_csv(f, header=param_header)
                with open(self.outpath_cycles, 'a') as f:
                    cycle_summary.to_csv(f, header=param_header)

        self.logger.info('Simulation completed without errors.')


'''
param_lst = [
    {'p_tandup': 0.000, 'DM_init': [1]*9+[2], 'resistance_fx':'linear'},
    {'p_tandup': 0.001, 'DM_init': [1]*9+[2], 'resistance_fx':'linear'},
    {'p_tandup': 0.010, 'DM_init': [1]*9+[2], 'resistance_fx':'linear'},
    {'p_tandup': 0.100, 'DM_init': [1]*9+[2], 'resistance_fx':'linear'},
    {'p_tandup': 1.000, 'DM_init': [1]*9+[2], 'resistance_fx':'linear'},
    {'p_tandup': 0.000, 'DM_init': [1], 'resistance_fx':'linear'},
    {'p_tandup': 0.001, 'DM_init': [1], 'resistance_fx':'linear'},
    {'p_tandup': 0.010, 'DM_init': [1], 'resistance_fx':'linear'},
    {'p_tandup': 0.100, 'DM_init': [1], 'resistance_fx':'linear'},
    {'p_tandup': 1.000, 'DM_init': [1], 'resistance_fx':'linear'},
    {'p_tandup': 0.000, 'DM_init': [1]*9+[2], 'resistance_fx':'asympt'},
    {'p_tandup': 0.001, 'DM_init': [1]*9+[2], 'resistance_fx':'asympt'},
    {'p_tandup': 0.010, 'DM_init': [1]*9+[2], 'resistance_fx':'asympt'},
    {'p_tandup': 0.100, 'DM_init': [1]*9+[2], 'resistance_fx':'asympt'},
    {'p_tandup': 1.000, 'DM_init': [1]*9+[2], 'resistance_fx':'asympt'},
    {'p_tandup': 0.000, 'DM_init': [1], 'resistance_fx':'asympt'},
    {'p_tandup': 0.001, 'DM_init': [1], 'resistance_fx':'asympt'},
    {'p_tandup': 0.010, 'DM_init': [1], 'resistance_fx':'asympt'},
    {'p_tandup': 0.100, 'DM_init': [1], 'resistance_fx':'asympt'},
    {'p_tandup': 1.000, 'DM_init': [1], 'resistance_fx':'asympt'}
]
param_lst = [
    {'p_tandup': 0.000, 'DM_init': [1]*7+[3], 'resistance_fx':'linear', 'strains': None, 'experiment': '[1]*7+[3]'},
    {'p_tandup': 0.000, 'DM_init': [1]*8+[2], 'resistance_fx':'linear', 'strains': None, 'experiment': '[1]*8+[2]'},
    {'p_tandup': 0.000, 'DM_init': [1]*10, 'resistance_fx':'linear', 'strains': None, 'experiment': '[1]*10'}
]
param_lst = [
    {'p_tandup': 0.000, 'DM_init': [1]*7+[3], 'resistance_fx':'linear', 'strains': None, 'experiment': '[1]*7+[3]'}
]

param_p_tandup_lst = [0, 0.001, 0.01, 0.1, 1.0]
sim_param_lst = []
dm_tbl_last_allparam = []
for sim_idx in range(0, len(param_lst)):
    print 'Running parameter setting {0} of {1}'.format(str(sim_idx+1), str(len(param_lst)))
    this_param = param_lst[sim_idx]
    sim_lst = []
    dm_tbl_lst = []
    for sim in range(0,num_simulations):
        print 'Running simulation number {0} of {1}'.format(str(sim+1), str(num_simulations))
        simulation = Simulation(num_cells=num_cells,
                                param_p_tandup=this_param['p_tandup'],
                                DM_init=this_param['DM_init'],
                                resistance_fx=this_param['resistance_fx'],
                                strains=this_param['strains'])
        simulation.run()
        cell_tbl, summary_tbl, dm_tbl = simulation.summarise()
        summary_tbl['dm_sums'].astype('int', inplace=True)
        summary_tbl['sim'] = sim
        sim_lst.append(summary_tbl)
        dm_tbl['sim'] = sim
        dm_tbl_lst.append(dm_tbl)
    sim_tbl = pd.concat(sim_lst)
    sim_tbl['experiment'] = this_param['experiment']
    sim_param_lst.append(sim_tbl)
    dm_tbl = pd.concat(dm_tbl_lst)
    dm_tbl['experiment'] = this_param['experiment']
    dm_tbl_last_allparam.append(dm_tbl)

pd.concat(sim_param_lst).to_csv('output/sim_tbl_multi.csv')
pd.concat(dm_tbl_last_allparam).to_csv('output/dm_tbl_multi.csv')


sim_counts = sim_tbl.groupby(['sim', 'dm_sums'], as_index=False).count().rename(columns={'id': 'counts'})
#sim_summary = pd.DataFrame({'id':pd.unique(sim_tbl['id'])})

#ax = sim_counts.boxplot(column='counts', by='dm_sums', return_type='axes')
#fig = ax.get_figure()
#fig.savefig('output/boxplot.png')
#plt.close(fig)

sim_tbl.to_csv('output/sim_tbl.csv')
'''
'''
fig = plt.figure(3)
plt.boxplot(sim_counts[['dm_sums', 'counts']])
#plt.plot(summary_tbl['dm_sums'], summary_tbl['dm_means'], 'k+')
plt.xlabel('Total cassettes')
plt.ylabel('Mean cassettes per DM')
fig.savefig('output/boxplot.png')
plt.close(fig)
'''

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

    return logger

# Tests
if False:
    c1 = Cell(DM_lst=[1,2,3])
    c2 = Cell()
    sum([c1,c2])

# Testrun
if __name__ == "__main__":
    sim_job = SimulationJob(path_json_params='input/params_exp2_debug.json',
                            outpath='output/exp2_debug/',
                            logger=generate_logger('output/exp2_debug/', 'test'),
                            job_id=None, random_seed=0)
    sim_job.run()
