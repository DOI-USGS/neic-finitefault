def parser_add_tensor(parser):
    """
    """
    parser.add_argument(
        "-gcmt", "--gcmt_tensor",
        help="location of GCMT moment tensor file")
    parser.add_argument(
        "-qcmt", "--qcmt_tensor",
        help="location of QuakeML moment tensor file")
    return parser


def parser_ffm_data(parser):
    """
    """
    parser.add_argument(
        "-t", "--tele", action="store_true",
        help="use teleseismic data in modelling")
    parser.add_argument(
        "-su", "--surface", action="store_true",
        help="use surface waves data in modelling")
    parser.add_argument(
        "-st", "--strong", action="store_true",
        help="use strong motion data in modelling")
    parser.add_argument(
        "--cgps", action="store_true",
        help="use cGPS data in modelling")
    parser.add_argument(
        "--gps", action="store_true",
        help="use GPS data in modelling")
    parser.add_argument(
        "-in", "--insar", action="store_true",
        help="use InSar data in modelling")
    return parser


def parser_data_process(parser):
    """
    """
    parser.add_argument(
        "-t", "--tele", action="store_true",
        help="process teleseismic data")
    parser.add_argument(
        "-su", "--surface", action="store_true",
        help="process surface waves data")
    parser.add_argument(
        "-st", "--strong", action="store_true",
        help="process strong motion data")
    parser.add_argument(
        "--cgps", action="store_true",
        help="process cGPS data")
    return parser


def parser_data_dict(parser):
    """
    """
    parser.add_argument(
        "-t", "--tele", action="store_true",
        help="create JSON for teleseismic body waves")
    parser.add_argument(
        "-su", "--surface", action="store_true",
        help="create JSON for surface waves")
    parser.add_argument(
        "-st", "--strong", action="store_true",
        help="create JSON for strong motion data")
    parser.add_argument(
        "--cgps", action="store_true",
        help="create JSON for cGPS data")
    parser.add_argument(
        "--gps", action="store_true",
        help="create JSON for static GPS data")
    parser.add_argument(
        "-in", "--insar", action="store_true",
        help="create JSON for InSar data")
    parser.add_argument(
        "-ina", "--insar_asc", default=[], nargs='*',
        help="Ascending InSar data")
    parser.add_argument(
        "-ind", "--insar_desc", default=[], nargs='*',
        help="Descending InSar data")
    parser.add_argument(
        "-inar", "--insar_asc_ramp", nargs='*',
        default=None, help="Ramp of ascending InSar data")
    parser.add_argument(
        "-indr", "--insar_desc_ramp", nargs='*',
        default=None, help="Ramp of descending InSar data")
    return parser


def parser_fill_data_files(parser):
    """
    """
    parser.add_argument(
        "-t", "--tele", action="store_true",
        help="fill text files with teleseismic data")
    parser.add_argument(
        "-su", "--surface", action="store_true",
        help="fill text files with surface waves data")
    parser.add_argument(
        "-st", "--strong", action="store_true",
        help="fill text files with strong motion data")
    parser.add_argument(
        "--cgps", action="store_true",
        help="fill text files with cGPS data")
    parser.add_argument(
        "--gps", action="store_true",
        help="fill text files with static GPS data")
    parser.add_argument(
        "-in", "--insar", action="store_true",
        help="fill text files with InSar data")
    return parser


def parser_data_plot(parser):
    """
    """
    parser.add_argument(
        "-t", "--tele", action="store_true",
        help="plot misfit of teleseismic data")
    parser.add_argument(
        "-su", "--surface", action="store_true",
        help="plot misfit of surface waves data")
    parser.add_argument(
        "-st", "--strong", action="store_true",
        help="plot strong motion stations and strong motion misfit")
    parser.add_argument(
        "--cgps", action="store_true",
        help="plot misfit of cGPS data")
    parser.add_argument(
        "--gps", action="store_true",
        help="plot GPS data")
    parser.add_argument(
        "-in", "--insar", action="store_true",
        help="plot InSar data")
    return parser


def parser_add_gf(parser):
    """
    """
    parser.add_argument(
        "-t", "--tele", action="store_true",
        help="compute teleseismic body waves GF")
    parser.add_argument(
        "-st", "--strong", action="store_true",
        help="compute strong motion GF")
    parser.add_argument(
        "--cgps", action="store_true",
        help="compute cGPS GF")
    parser.add_argument(
        "--gps", action="store_true",
        help="compute static GPS GF")
    parser.add_argument(
        "-in", "--insar", action="store_true",
        help="get static INSAR GF")
    return parser


def get_used_data(args):
    """
    """
    used_data = []
    if 'gps' in args:
        used_data = used_data + ['gps'] if args.gps else used_data
    if 'strong' in args:
        used_data = used_data + ['strong_motion'] if args.strong else used_data
    if 'cgps' in args:
        used_data = used_data + ['cgps'] if args.cgps else used_data
    if 'tele' in args:
        used_data = used_data + ['tele_body'] if args.tele else used_data
    if 'surface' in args:
        used_data = used_data + ['surf_tele'] if args.surface else used_data
    if 'insar' in args:
        used_data = used_data + ['insar'] if args.insar else used_data
    return used_data
