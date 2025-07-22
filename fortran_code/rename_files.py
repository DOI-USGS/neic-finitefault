import os
import glob
from shutil import move


REPLACE_DICT = {
    'Readlp.das': 'channels_body.txt',
    'Obser.tele': 'waveforms_body.txt',
    'Wave.tele': 'wavelets_body.txt',
    'Weight': 'body_wave_weight.txt',
    'instrumental_response': 'instrumental_response.txt',
    'synm.tele': 'synthetics_body.txt',
    'Readlp.inf_low': 'channels_surf.txt',
    'Obser.str_low': 'waveforms_surf.txt',
    'Wave.str_low': 'wavelets_surf.txt',
    'synm.str_low': 'synthetics_surf.txt',
    'Readlp.inf': 'channels_strong.txt',
    'Obser.str': 'waveforms_strong.txt',
    'Wave.str': 'wavelets_strong.txt',
    'synm.str': 'synthetics_strong.txt',
    'Readlp.cgps': 'channels_cgps.txt',
    'Obser.cgps': 'waveforms_cgps.txt',
    'Wave.cgps': 'wavelets_cgps.txt',
    'synm.cgps': 'synthetics_cgps.txt',
    'filter_tele': 'filter_tele.txt',
    'filter_strong': 'filter_strong.txt',
    'Readlp.static': 'static_data.txt',
    'synm.static': 'static_synthetics.txt',
    'Fault.time': 'fault&rise_time.txt',
    'Fault.pos': 'point_sources.txt',
    'Niu_model': 'shear_model.txt',
    'vel_model': 'vel_model.txt',
    'Green_static_subfault': 'Green_static_subfault.txt',
    'HEAT.IN': 'annealing.txt',
    'bound.in': 'model_space.txt',
    'bound.special': 'special_model_space.txt',
    'continue': 'regularization_borders.txt',
    'continue.special': 'special_regularization_borders.txt',
    'Green.in': 'Green_strong.txt',
    'Green_cgps.in': 'Green_cgps.txt'
}


def replace():
    """
    """
    files = os.listdir()
    files = [v for v in files if os.path.isfile(v)]
    for file in files:
        if file in REPLACE_DICT:
            move(file, REPLACE_DICT[file])


if __name__ == '__main__':
    replace()
