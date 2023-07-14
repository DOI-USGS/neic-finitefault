import argparse

from wasp.manage_parser import (
    get_used_data,
    parser_add_gf,
    parser_add_tensor,
    parser_data_dict,
    parser_data_plot,
    parser_data_process,
    parser_ffm_data,
    parser_fill_data_files,
)


def test_load_parsers():
    parser_add_tensor(argparse.ArgumentParser())
    parser_ffm_data(argparse.ArgumentParser())
    parser_data_process(argparse.ArgumentParser())
    parser_data_dict(argparse.ArgumentParser())
    parser_fill_data_files(argparse.ArgumentParser())
    parser_data_plot(argparse.ArgumentParser())
    parser_add_gf(argparse.ArgumentParser())


def test_get_used_data():
    args = argparse.Namespace(gps=True, tele=True, insar=True)
    assert get_used_data(args) == ["gps", "tele_body", "insar"]
    args = argparse.Namespace(strong=True, cgps=True, surface=True)
    assert get_used_data(args) == ["strong_motion", "cgps", "surf_tele"]
