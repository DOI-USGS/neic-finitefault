import os
import json
import glob
import management as mng
import logging
import management as mng
from obspy.io.sac import SACTrace


def get_sacpz_file(sacheader, data='obspy'):
    """
    """
    loc = sacheader.khole
    if data == 'IRIS' and not loc:
        loc = '--'
    else:
        if not loc:
            loc = '__'
        elif len(loc) == 0:
            loc = '__'
    pzfile = 'SAC_PZs_{}_{}_{}_{}'.format(sacheader.knetwk, sacheader.kstnm,
                                          sacheader.kcmpnm, loc)
    if data == 'IRIS':
        pzfile = 'SACPZ.{}.{}.{}.{}'.format(sacheader.knetwk, sacheader.kstnm,
                                            loc, sacheader.kcmpnm)
    return pzfile


def convert_response_acc(resp_file):
    """
    """
    with open(resp_file, 'r') as infile:
        lines = [line for line in infile]
    indexes = [i for i, line in enumerate(lines) if 'ZEROS' in line.split()]
    index0 = 0
    lines2 = []
    for index in indexes:
        lines2 = lines2 + lines[index0:index]
        string, nzeroes = lines[index].split()
        if nzeroes == '2':
            lines2 = lines2 + ['ZEROS\t 0\n']
            index0 = index + 3
        if nzeroes == '3':
            lines2 = lines2 + ['ZEROS\t 0\n']
            index0 = index + 4
        lines2 = lines2 + lines[index0:]
    with open(resp_file, 'w')as outfile:
        for line in lines2:
            outfile.write(line)
    return resp_file


def replace_response_sac(sac_files, old_responses, names, new_response=None,
                         pre_processing=None, add_response=False,
                         freq0=None, freq1=None, freq2=None, freq3=None):
    """
    """
    zipped = zip(sac_files, old_responses, names)
    input_sac = ''
    if not pre_processing:
        pre_processing = 'rmean \n'
    bandpass = ''
    if freq0 and freq1 and freq2 and freq3:
        bandpass = ' freq {} {} {} {}\n'.format(freq0, freq1, freq2, freq3)
    if add_response:
        post_processing = ' to polezero s {}{}\n'.format(new_response, bandpass)
    else:
        post_processing = '{}\n'.format(bandpass)
    instructions = lambda resp0:\
        pre_processing + 'transfer from polezero s {}'.format(resp0)\
        + post_processing
    for sac_file, old_resp, name in zipped:
        input_sac = '{}\nread {} \n'.format(input_sac, sac_file)\
                    + instructions(old_resp)\
                    + 'write {}'.format(name)

    input_sac = '{}\nquit\n'.format(input_sac)
    out, err = mng.run_sac(input_sac)
    return out, err


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    args = parser.parse_args()
    os.chdir(args.folder)
    tensor_info = {
        'depth': 25,
        'time_shift': 40
    }
    data_prop = {
        'tele_filter': {
            'freq0': 0.001,
            'freq1': 0.002,
            'freq2': 1.0,
            'freq3': 1.2
        }
    }
    files0 = glob.glob('SACPZ*')
    for file in files0:
        convert_response_acc(file)
    files = glob.glob('*SAC')
    __remove_response_str2(files, data='IRIS')