#!/usr/bin/env python3

SUPPLY_VOLTAGE   = 12.19
SHUNT_RESISTANCE = 0.02
TICKS_PER_SECOND = 1735786190
SAMPLES_PER_SECOND = 10
IDLE_AMPERAGE = 0.004
CLOCK_VARIANCE_CUTOFF = 2.1
PAD = 100

IDLE_HEAD = 14
IDLE_TAIL = 14


import numpy
import math
import subprocess as SP
import sys
import glob
import os.path as path
from pprint import pprint


def func_from_filename(filename):
    filename = path.basename(filename)
    filename = filename.upper()
    parts = filename.split('_')
    if parts[0] != "IAN": # REMOVE BEFORE PUTTING IN SVN
        print(filename)
    assert(parts[-1][-3:] == "CSV")
    assert(parts[-2] == "RUN")
    parts = parts[1:-2]
    mod = parts[-1]
    func = '_'.join(parts[:-1])
    return func


def split_func(func):
    parts = func.split('_')
    if parts[-1] in ["MIX", "ADAPTIVE"]:
        typ = parts[-1]
        func = '_'.join(parts[:-1])
    else:
        typ = '_'.join(parts[-2:])
        func = '_'.join(parts[:-2])
    return (func, typ)


def msum(iterable):
    "Full precision summation using multiple floats for intermediate values"
    # Rounded x+y stored in hi with the round-off stored in lo.  Together
    # hi+lo are exactly equal to x+y.  The inner loop applies hi/lo summation
    # to each partial so that the list of partial sums remains exact.
    # Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:
    # www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
    partials = []               # sorted, non-overlapping partial sums
    for x in iterable:
        i = 0
        for y in partials:
            if abs(x) < abs(y):
                x, y = y, x
            hi = x + y
            lo = y - (hi - x)
            if lo:
                partials[i] = lo
                i += 1
            x = hi
        partials[i:] = [x]
    return sum(partials, 0.0)


def read_amperage(filename):
    with open(filename, 'r') as f:
        data = [line.split() for line in f.readlines() if line.strip() != ""]
    data = [(int(d[0]), float(d[1])) for d in data]
    return data


def amperage_to_power(amperage):
    return [I*SUPPLY_VOLTAGE for I in amperage]


def power_to_energy(power, timestamps):
    energy = list()
    for i in range(0, len(power)-1):
        delta_ticks = timestamps[i+1] - timestamps[i]
        delta_time = delta_ticks / TICKS_PER_SECOND
        assert(delta_time > 0)
        e = delta_time * (power[i]+power[i+1])/2
        energy.append(e)

    return energy


def check_variance(timestamps):
    # how long between each sample
    time_diffs = [timestamps[i]-timestamps[i-1] for i in range(1, len(timestamps))]
    avg_sample_length = msum(time_diffs)/len(time_diffs)

    # the diffence between adjacent sample lengths 
    time_diff_diffs = [abs(time_diffs[i]-time_diffs[i-1]) for i in range(1, len(time_diffs))]

    max_var_between_samples = (max(time_diff_diffs)+avg_sample_length)/avg_sample_length
    return max_var_between_samples


def check_idle(amperage):
    # was an OS operation going on at the start?
    if msum(amperage[0:IDLE_HEAD]) / IDLE_HEAD > IDLE_AMPERAGE:
        return True
    # was an OS operation going on at the end?
    return msum(amperage[0:IDLE_TAIL]) / IDLE_TAIL > IDLE_AMPERAGE


def align_start_times(function_dict):
    for func in function_dict.keys():
        for typ in function_dict[func].keys():
            for i in range(len(function_dict[func][typ])):
                # remove idle cost
                this_idle = msum(function_dict[func][typ][i][:IDLE_HEAD])/IDLE_HEAD
                function_dict[func][typ][i] = [d-this_idle for d in function_dict[func][typ][i]]

                # set upswing to be the IDLE_HEADth item
                this_idle_max = max(function_dict[func][typ][i][:IDLE_HEAD])
                for j in range(len(function_dict[func][typ][i])):
                    if function_dict[func][typ][i][j] > this_idle_max*3:
                        if j - IDLE_HEAD < 3:
                            break
                        start = j - IDLE_HEAD
                        assert(start >= 0)
                        function_dict[func][typ][i] = function_dict[func][typ][i][start:]
                        break

                for j in range(2*IDLE_HEAD, len(function_dict[func][typ][i])):
                    if function_dict[func][typ][i][j] < this_idle_max*3:
                        end = j + IDLE_TAIL
                        function_dict[func][typ][i] = function_dict[func][typ][i][:end]
                        break
                       

def fft_pass(amperage):
    arr = [0 for _ in range(PAD)] + amperage + [0 for _ in range(PAD)]

    res = numpy.fft.fft(arr)

    mags = []

    for r in res:
        mags.append(2.0/len(arr)*abs(r))

    freqs = []
    for i in range(len(arr)):
        freqs.append(i*SAMPLES_PER_SECOND/len(arr))

    M = max(mags)
    base_name = sys.argv[1][:-4]

    for i in range(len(mags)):
        if freqs[i] > 40 and freqs[i] < 130:
            res[i] = complex(0)

                    
    reverse = numpy.fft.ifft(res)
    return [abs(I) for I in reverse[PAD:-PAD]]

    # with open(base_name + "_freq.csv", 'w') as f:
    #     for i in range(min(len(freqs), len(mags), len(arr))):
    #         f.write("{},{}\n".format(freqs[i], mags[i]))

                   #res[i] = complex(res[i])*complex(mags[i]/M)
                
                

                
def read_power_directory(dirname):
    energy_dict = dict()
    power_dict = dict()
    data = list()
    for filename in sorted(glob.glob(dirname+"/*.CSV")):
        func = func_from_filename(filename)
        func, typ = split_func(func)
        if func not in energy_dict:
            energy_dict[func] = dict()
            power_dict[func] = dict()
        if typ not in energy_dict[func]:
            energy_dict[func][typ] = list()
            power_dict[func][typ] = list()
            
        new_datum = read_amperage(filename)
        if len(new_datum) == 0:
            print("BAD FILE: {}".format(filename))
            continue
        amperage = [d[1] for d in new_datum]
        amperage = fft_pass(amperage)
        timestamps = [d[0] for d in new_datum]
        var = check_variance(timestamps)
        if var > CLOCK_VARIANCE_CUTOFF:
            print("Dropping data due to clock varance of: {}".format(var))
            print("Datafile: {}".format(filename))
            continue
        if check_idle(amperage):
            print("Dropping run due to high start amperage")
            print("Datafile: {}".format(filename))
            continue
        power = amperage_to_power(amperage)
        power_dict[func][typ].append(power)
        energy = power_to_energy(power, timestamps)
        energy_dict[func][typ].append(energy)


    align_start_times(energy_dict)
    align_start_times(power_dict)
    return energy_dict, power_dict


def read_logfile(filename):
    with open(filename) as f:
        lines = f.readlines()
    time_dict = dict()
    iter_dict = dict()
    func = None
    typ = None
    for i in range(len(lines)):
        if lines[i].strip() == "":
            continue
        parts = lines[i].replace(" ","").split(",")
        if parts[0] == "Iterations":
            iters = int(parts[1])
            func = func_from_filename(lines[i+1].replace(" ","").split(',')[0])
            func, typ = split_func(func)
            if func not in iter_dict:
                iter_dict[func] = dict()
                time_dict[func] = dict()
            if typ not in iter_dict[func]:
                iter_dict[func][typ] = iters
                time_dict[func][typ] = list()
            else:
                print("Repeated data in logfile for fun '{}' type '{}'".format(func, typ))
                sys.exit(-1)
            continue
        
        func = func_from_filename(parts[0])
        func, typ = split_func(func)
        time_dict[func][typ].append(int(parts[1]))

    return time_dict, iter_dict


def dump_energy(energy_matrix, filename, headerp):
    with open(filename+".dat", 'w') as f:
        f.write("runs\t")
        if not headerp:
            for i in range(len(energy_matrix)):
                f.write("run_{}\t".format(i))
            f.write("\n")
            
        for i in range(max([len(col) for col in energy_matrix])):
            f.write("{}\t".format(i/SAMPLES_PER_SECOND))
            for col in energy_matrix:
                try:
                    f.write("{}\t".format(col[i]))
                except:
                    f.write("#\t")
            f.write("\n")
                    

def config_plot(xmax, ymax, filename, plotname, keyp):
    with open(filename+".gnu", 'w') as f:
        f.write("set output '{}.png'\n".format(plotname))
        f.write("set title '{}'\n".format(plotname))
        f.write("set term png size 1920,1080\n")
        f.write("set datafile missing '#'\n")
        f.write("set style data lines\n")
        if not keyp:
            f.write("set key off\n")
        f.write("set xlabel 'Time (Seconds)'\n")
        f.write("set ylabel 'Power Above Idle (Watts)'\n")
        f.write("set yrange [-0.01:{}]\n".format(ymax))
        f.write("set xrange [0:{}]\n".format(xmax))
        f.write("file = '{}.dat\n".format(filename))
        f.write("cols = int(system('head -1 '.file.' | wc -w'))\n")
        f.write("plot for [i=2:cols] '{}.dat' using 1:i title columnheader(i+1)\n".format(filename))
        

def gnuplot(filename):
    command = "gnuplot "+filename+".gnu"
    with SP.Popen(command,
                  stdout=SP.PIPE, stderr=SP.STDOUT,shell=True) as proc:
        output = proc.stdout.read().decode("utf-8")
        proc.wait()

        if proc.returncode != 0:
            print("Unable to run gnuplot")
            print("Non-zero return code: {}".format(proc.returncode))
            print("Command used: {}".format(command))
            print("Trace:\n{}".format(output))
            sys.exit(proc.returncode)

    return output


def average_energy(energy):
    avg = list()
    for i in range(max([len(col) for col in energy])):
        e = list()
        for col in energy:
            try:
                e.append(col[i])
            except:
                pass
        if len(e) < 0.5*len(energy):
            break
        avg.append(msum(e)/len(e))
        #avg.append(sorted(e)[len(e)//2])
    return avg


def plot(func_name, type_dict):
    print("plotting {}".format(func_name))
    xmax_list = [max([len(col) for col in em]) for em in type_dict.values()]
    xmax = math.ceil(max(xmax_list)/10)*10/SAMPLES_PER_SECOND

    ymax_list = [max([max(col) for col in em]) for em in type_dict.values()]
    ymax = math.ceil(max(ymax_list)*1000) / 1000

    comparative_data = list()
    for typ in sorted(type_dict.keys()):
        avg_energy = average_energy(type_dict[typ])
        comparative_data.append([typ] + fft_pass(avg_energy))
        #dump_energy(type_dict[typ], "temp", False)
        #config_plot(xmax, ymax, "temp", "{}_{}".format(func_name, typ), False)
        #gnuplot("temp")


    dump_energy(comparative_data, "temp", True)
    config_plot(xmax, ymax, "temp", "{}".format(func_name), True)
    gnuplot("temp")
    
    print("done\n")

        
def main():
    if len(sys.argv) != 3:
        print("usage: {} logfile directory".format(sys.argv[0]))
        return
    energy_dict, amperage_dict = read_power_directory(sys.argv[1])
    time_dict, iterations_dict = read_logfile(sys.argv[2])

    sep = ","
    print(sep.join(["Function", "Precision", "Energy Per Iteration (Joules)", "Time Per Iteration (Milliseconds)", "Energy_Relative To 64 Bit", "Time Relative To 64 Bit"]))
    for func in sorted(energy_dict.keys()):
        lines = list()
        energy64 = None
        time64 = None
        for typ in energy_dict[func].keys():
            iters = iterations_dict[func][typ]

            total_energy = msum([msum(col) for col in energy_dict[func][typ]])
            energy_per_run = total_energy / len(energy_dict[func][typ])
            energy_per_iter = energy_per_run / iters
            
            total_time = sum(time_dict[func][typ])
            time_per_run = total_time / len(time_dict[func][typ])
            time_per_iter = time_per_run / iters
            
            lines.append([func, typ, str(energy_per_iter), str(time_per_iter)])
            if "64" in typ:
                energy64 = energy_per_iter
                time64 = time_per_iter
        for line in lines:
            line.append(str(float(line[2])/energy64))
            line.append(str(float(line[3])/time64))
            print(sep.join(line))
            
        print("")

    # Do plotting
    for func in sorted(amperage_dict.keys()):
        plot(func, amperage_dict[func])
            
    
if __name__ == "__main__":
    main()
