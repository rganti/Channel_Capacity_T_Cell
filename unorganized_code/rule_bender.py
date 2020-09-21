import os

from realistic_network import TcrCycleSelfWithForeign, TcrCycleSelfLigand
from two_species import KPSingleSpecies


class RBTcrSelfForeign(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)

    def modify_forward_reverse(self, reactants, products, forward_rate, reverse_rate):
        forward_key = ''.join(reactants) + '_' + ''.join(products)
        self.forward_rates[forward_key] = forward_rate
        self.forward_rxns.append([reactants, products])

        reverse_key = ''.join(products) + '_' + ''.join(reactants)
        self.reverse_rates[reverse_key] = reverse_rate
        self.reverse_rxns.append([products, reactants])

    def create_loop(self, reactants, products, rate):
        forward_key = ''.join(reactants) + '_' + ''.join(products)
        self.forward_rates[forward_key] = rate
        self.forward_rxns.append([reactants, products])


class KPParameterTest(KPSingleSpecies):
    def __init__(self, self_foreign=False, arguments=None):
        KPSingleSpecies.__init__(self, self_foreign=self_foreign, arguments=arguments)

        if self.self_foreign_flag:
            print("initializing Self and Foreign")
            self.ligand = RBTcrSelfForeign(arguments=arguments)
        else:
            print("initializing Self Only")
            self.ligand = TcrCycleSelfLigand(arguments=arguments)

        self.ligand.n_initial["Ls"] = 800
        if self.arguments.run:
            self.run_time = self.arguments.run
        else:
            self.run_time = 1000

        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

        self.simulation_time = 20


class RuleBenderSharedCommands(object):

    def __init__(self, n_initial, record, script_name):
        self.n_initial = n_initial
        self.record = record
        self.script_name = script_name

    def initialize_species(self, f):
        f.write("begin seed species\n")
        for key, value in self.n_initial.items():
            f.write("\t{0} {1}\n".format(key, value))
        f.write("end seed species\n\n")

    def record_species(self, f):
        f.write("begin observables\n")
        for item in self.record:
            f.write("Molecules {0} {0}()\n".format(item))
        f.write("end observables\n\n")


class KPRuleBender(object):
    def __init__(self, self_foreign=False, arguments=None):
        self.self_foreign = self_foreign
        self.arguments = arguments

        if self.self_foreign:
            print("initializing Self and Foreign")
            self.ligand = TcrCycleSelfWithForeign(arguments=arguments)
        else:
            print("initializing Self Only")
            self.ligand = TcrCycleSelfLigand(arguments=arguments)

        self.ligand.n_initial["Ls"] = 800

        self.num_files = 100

        if self.arguments.run:
            self.run_time = self.arguments.run
        else:
            self.run_time = 1000

        self.simulation_time = 10
        self.home_directory = os.getcwd()

        self.num_kp_steps = 1
        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

    def set_time_step(self):
        if self.arguments:
            if self.arguments.ss:
                time_step = 1.0
            else:
                time_step = self.run_time
        else:
            time_step = self.run_time
        return time_step

    @staticmethod
    def define_reactions(f, rxn, rate):
        for i in range(len(rxn)):
            input_string = "rxn "
            destroy_string = ""
            input = rxn[i][0]
            output = rxn[i][1]

            rate_key = ""
            for item in input:
                input_string += "{0}:{1} ".format(item.lower(), item)
                destroy_string += "destroy {0}; ".format(item.lower())
                rate_key += item

            rate_key += "_"
            output_string = ""
            for item in output:
                output_string += "new {0}; ".format(item)
                rate_key += item

            rate_string = "at {0} -> ".format(rate[rate_key])

            f.write(input_string + rate_string + destroy_string + output_string + "\n")

    def generate_rule_bender_script(self, simulation_name):
        script_name = simulation_name + ".bngl"
        shared = RuleBenderSharedCommands(self.ligand.n_initial, self.record, script_name)

        f = open(script_name, "w")

        shared.initialize_species(f)
        shared.record_species(f)

        self.define_reactions(f, self.forward_rxns, self.forward_rates)
        self.define_reactions(f, self.reverse_rxns, self.reverse_rates)

        f.close()

    def generate_qsub(self, simulation_name, time_step):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime={1},nodes=1:ppn=1 -N {0}\n\n".format(simulation_name,
                                                                       datetime.timedelta(
                                                                           minutes=self.simulation_time)))
        q.write("cd $PBS_O_WORKDIR\n\n")
        q.write("echo $PBS_JOBID > job_id\n")
        q.write("EXE_FILE={0}\n".format(simulation_name))
        q.write("RUN_TIME={0}\n".format(self.run_time))
        q.write("STEP={0}\n\n".format(time_step))
        q.write("for j in {1.." + str(self.num_files) + "}\n")
        q.write("do\n")
        if time_step == self.run_time:
            q.write("\t ./$EXE_FILE -e $RUN_TIME > traj_$j\n")
        else:
            q.write("\t ./$EXE_FILE -e $RUN_TIME -t $STEP > traj_$j\n")
        q.write("done\n\n")
        q.write("python ~/SSC_python_modules/post_process.py --num_files {0} "
                "--run_time {1} --time_step {2}\n".format(self.num_files, self.run_time, time_step))
        if self.arguments.ss:
            q.write("python ~/SSC_python_modules/plot.py \n")
        q.close()

    def generate(self, simulation_name, time_step):
        self.generate_ssc_script(simulation_name)
        compile_script(simulation_name + ".rxn")
        self.generate_qsub(simulation_name, time_step)

    def single_add_step(self):
        self.num_kp_steps += 1
        self.ligand.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))

        self.forward_rates[self.ligand.symbol + "{0}".format(self.num_kp_steps - 2)] = self.ligand.rate_constants.k_p
        self.forward_rxns.append([[self.ligand.symbol + "{0}".format(self.num_kp_steps - 2)],
                                  [self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)]])

        self.reverse_rates[self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)] = self.reverse_rates[
            self.ligand.symbol + "0"]
        self.reverse_rxns.append([[self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)], self.ligand.inputs])

        self.record.append(self.ligand.symbol + "{0}".format(self.num_kp_steps - 1))

    def competing_add_step(self):
        self.num_kp_steps += 1
        self.ligand.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))
        for i in ["C", "D"]:
            self.forward_rates[i + "{0}".format(self.num_kp_steps - 2)] = self.ligand.rate_constants.k_p
            self.forward_rxns.append([[i + "{0}".format(self.num_kp_steps - 2)],
                                      [i + "{0}".format(self.num_kp_steps - 1)]])
            if i == "C":
                self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.ligand.rate_constants.k_foreign_off
                self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Lf"]])
            elif i == "D":
                self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.ligand.rate_constants.k_self_off
                self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Ls"]])
            self.record.append(i + "{0}".format(self.num_kp_steps - 1))

    def add_step(self):
        if self.self_foreign_flag:
            self.competing_add_step()
        else:
            self.single_add_step()

    def main_script(self, run=False):
        sample = []
        for i in range(self.ligand.num_samples):
            directory = "sample_" + str(i)
            s = self.ligand.p_ligand[i]
            sample.append(s)
            self.ligand.change_ligand_concentration(s)
            simulation_name = self.ligand.simulation_name + "_" + str(i)

            os.makedirs(directory)
            print("Made " + directory)
            os.chdir(directory)
            print("Changed into directory: " + str(os.getcwd()))

            if self.ligand.num_kp_steps > 6:
                self.run_time = 1000

            self.generate(simulation_name, self.set_time_step())
            if run:
                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()
            os.chdir(self.home_directory)

        np.savetxt("Ligand_concentrations", sample, fmt='%f')
        np.savetxt("Ligand_concentrations_sorted", np.sort(sample), fmt='%f')
