import pathlib
from shutil import rmtree
from tempfile import mkdtemp
from ffm.model.result.solution import Solution
import unittest

SINGLE_SEGMENTS_DATA = {
    "rise_time": {"delta_rise": 1.5, "min_rise": 1.5, "windows": 10},
    "segments": [
        {
            "delta_dip": 14,
            "delta_strike": 17,
            "dip": 19,
            "dip_subfaults": 9,
            "hyp_dip": 5,
            "hyp_stk": 9,
            "neighbours": [],
            "rake": 109,
            "rupture_vel": 2.5,
            "stk_subfaults": 23,
            "strike": 6,
            "max_vel": 3.125,
            "min_vel": 1,
        }
    ],
}
MULTIPLE_SEGMENTS_DATA = {
    "rise_time": {"delta_rise": 1.5, "min_rise": 1.5, "windows": 10},
    "segments": [
        {
            "delay_segment": 0,
            "delta_dip": 14,
            "delta_strike": 17,
            "dip": 19,
            "dip_subfaults": 9,
            "hyp_dip": 5,
            "hyp_stk": 9,
            "neighbours": [],
            "rake": 109,
            "rupture_vel": 2.5,
            "stk_subfaults": 23,
            "strike": 6,
            "max_vel": 3.125,
            "min_vel": 1,
        },
        {
            "delay_segment": 1,
            "delta_dip": 12,
            "delta_strike": 18,
            "dip": 22,
            "dip_subfaults": 1,
            "hyp_dip": 5,
            "hyp_stk": 7,
            "neighbours": [],
            "rake": 105,
            "rupture_vel": 2.2,
            "stk_subfaults": 21,
            "strike": 7,
            "max_vel": 3.2,
            "min_vel": 2,
        },
    ],
    "connections": [
        {
            "segment1": 1,
            "segment2": 2,
            "segment1_subfault": [2, 1],
            "segment2_subfault": [2, 3],
        }
    ],
}

SINGLE_SEGMENT = """ #Total number of fault_segments=           1
#Fault_segment =   1 nx(Along-strike)= 23 Dx = 17.93km ny(downdip)=   9 Dy = 14.92km
#Boundary of Fault_segment    1. EQ in cell ( 9, 5). Lon:  -71.6700   Lat:  -31.5700
 #Lon.  Lat.  Depth
  -72.5035553      -32.8519516      0.716799974    
  -71.2028580      -32.9804459       44.0831985    
  -70.7055817      -29.3264771       44.0831985    
  -72.0062790      -29.1979809      0.716799974    
  -72.5035553      -32.8519516      0.716799974    
 #Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo
    -32.851952    -72.503555      0.716800    126.725410     89.414207      6.613912     19.280828     84.339615     13.500000      1.500000  0.119633E+27
    -32.691689    -72.481743      0.716800    144.565353     91.090263      6.613912     19.280828     75.787537      9.000000      1.500000  0.136475E+27
"""

MULTIPLE_SEGMENTS = """ #Total number of fault_segments=           2
#Fault_segment =   1 nx(Along-strike)=  7 Dx = 10.00km ny(downdip)=   5 Dy =  8.00km
#Boundary of Fault_segment    1. EQ in cell ( 3, 3). Lon:  -67.6280   Lat:   10.5460
 #Lon.  Lat.  Depth
  -68.7431793       10.4013395       2.26938462    
  -68.7114868       10.3201580       38.3306160    
  -68.1419678       10.5351620       38.3306160    
  -68.1736526       10.6163425       2.26938462    
  -68.7431793       10.4013395       2.26938462    
 #Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo
     10.401340    -68.743179      2.269385     76.072807    165.796600     69.000000     75.000000     12.644999      9.000000      3.000000  0.214500E+26
     10.433590    -68.657745      2.269385     85.640694    165.888412     69.000000     75.000000      8.947185     10.500000      4.500000  0.241478E+26
#Fault_segment =   2 nx(Along-strike)= 20 Dx = 11.00km ny(downdip)=   4 Dy =  7.00km
#Boundary of Fault_segment    2. EQ in cell ( 3, 3). Lon:  -67.6280   Lat:   10.5460
 #Lon.  Lat.  Depth
  -67.8415298       10.5379143       7.42419767    
  -67.8355103       10.4511595       43.4854279    
  -66.0394135       10.5715370       43.4854279    
  -66.0454254       10.6582918       7.42419767    
  -67.8415298       10.5379143       7.42419767    
 #Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo
     10.537914    -67.841530      7.424198     43.059509    200.494431     86.099998     75.000000     32.200001     12.000000      9.000000  0.121413E+26
     10.544036    -67.750198      7.424198    108.274605    183.671906     86.099998     75.000000     27.400000      1.500000      4.500000  0.305298E+26

"""


class TestSolution(unittest.TestCase):
    # runs before each test
    def setUp(self):
        self.tempdir = pathlib.Path(mkdtemp())
        self.solution_filepath = self.tempdir / "Solution.txt"

    # runs after each test even if test fails
    def tearDown(self):
        rmtree(self.tempdir)

    def test_parse_static_solution_multi_segment(self):
        # write sample single segment
        with open(self.solution_filepath, "w") as f:
            f.write(MULTIPLE_SEGMENTS)

        # parse solution
        solution = Solution.from_static_solution_file(
            segments_data=MULTIPLE_SEGMENTS_DATA,
            solution_filepath=self.solution_filepath,
        )

        # check values
        self.assertEqual(solution.num_segments, 2)
        self.assertEqual(solution.num_segments, len(solution.segments))
        # segment 1
        self.assertDictEqual(
            solution.connections[0].model_dump(),
            {
                "segment1": 1,
                "segment2": 2,
                "segment1_subfault": [2, 1],
                "segment2_subfault": [2, 3],
            },
        )
        segment1 = solution.segments[0]
        self.assertListEqual(list(segment1.data.lat), [10.401340, 10.433590])
        self.assertListEqual(list(segment1.data.lon), [-68.743179, -68.657745])
        self.assertListEqual(list(segment1.data.depth), [2.269385, 2.269385])
        self.assertListEqual(list(segment1.data.slip), [76.072807, 85.640694])
        self.assertListEqual(list(segment1.data.rake), [165.796600, 165.888412])
        self.assertListEqual(list(segment1.data.strike), [69, 69])
        self.assertListEqual(list(segment1.data.dip), [75, 75])
        self.assertListEqual(list(segment1.data.trup), [12.644999, 8.947185])
        self.assertListEqual(list(segment1.data.trise), [9, 10.5])
        self.assertListEqual(list(segment1.data.tfall), [3, 4.5])
        self.assertListEqual(list(segment1.data.moment), [0.214500e26, 0.241478e26])
        self.assertListEqual(
            segment1.polygon,
            [
                (-68.7431793, 10.4013395, 2.26938462),
                (-68.7114868, 10.3201580, 38.3306160),
                (-68.1419678, 10.5351620, 38.3306160),
                (-68.1736526, 10.6163425, 2.26938462),
                (-68.7431793, 10.4013395, 2.26938462),
            ],
        )
        self.assertEqual(segment1.nx, 7)
        self.assertEqual(segment1.dx, 10)
        self.assertEqual(segment1.ny, 5)
        self.assertEqual(segment1.dy, 8)
        metadata = segment1.metadata
        self.assertEqual(metadata.delay_segment, 0)
        self.assertEqual(metadata.delta_dip, 14)
        self.assertEqual(metadata.delta_strike, 17)
        self.assertEqual(metadata.dip, 19)
        self.assertEqual(metadata.dip_subfaults, 9)
        self.assertEqual(metadata.hyp_dip, 5)
        self.assertEqual(metadata.hyp_stk, 9)
        self.assertEqual(metadata.neighbours, [])
        self.assertEqual(metadata.rake, 109)
        self.assertEqual(metadata.rupture_vel, 2.5)
        self.assertEqual(metadata.stk_subfaults, 23)
        self.assertEqual(metadata.strike, 6)
        self.assertEqual(metadata.max_vel, 3.125)
        self.assertEqual(metadata.min_vel, 1)
        # segment 2
        segment2 = solution.segments[1]
        self.assertListEqual(list(segment2.data.lat), [10.537914, 10.544036])
        self.assertListEqual(list(segment2.data.lon), [-67.841530, -67.750198])
        self.assertListEqual(list(segment2.data.depth), [7.424198, 7.424198])
        self.assertListEqual(list(segment2.data.slip), [43.059509, 108.274605])
        self.assertListEqual(list(segment2.data.rake), [200.494431, 183.671906])
        self.assertListEqual(list(segment2.data.strike), [86.099998, 86.099998])
        self.assertListEqual(list(segment2.data.dip), [75, 75])
        self.assertListEqual(list(segment2.data.trup), [32.200001, 27.4])
        self.assertListEqual(list(segment2.data.trise), [12, 1.5])
        self.assertListEqual(list(segment2.data.tfall), [9, 4.5])
        self.assertListEqual(list(segment2.data.moment), [0.121413e26, 0.305298e26])
        self.assertListEqual(
            segment2.polygon,
            [
                (-67.8415298, 10.5379143, 7.42419767),
                (-67.8355103, 10.4511595, 43.4854279),
                (-66.0394135, 10.5715370, 43.4854279),
                (-66.0454254, 10.6582918, 7.42419767),
                (-67.8415298, 10.5379143, 7.42419767),
            ],
        )
        self.assertEqual(segment2.nx, 20)
        self.assertEqual(segment2.dx, 11)
        self.assertEqual(segment2.ny, 4)
        self.assertEqual(segment2.dy, 7)
        metadata2 = segment2.metadata
        self.assertEqual(metadata2.delay_segment, 1)
        self.assertEqual(metadata2.delta_dip, 12)
        self.assertEqual(metadata2.delta_strike, 18)
        self.assertEqual(metadata2.dip, 22)
        self.assertEqual(metadata2.dip_subfaults, 1)
        self.assertEqual(metadata2.hyp_dip, 5)
        self.assertEqual(metadata2.hyp_stk, 7)
        self.assertEqual(metadata2.neighbours, [])
        self.assertEqual(metadata2.rake, 105)
        self.assertEqual(metadata2.rupture_vel, 2.2)
        self.assertEqual(metadata2.stk_subfaults, 21)
        self.assertEqual(metadata2.strike, 7)
        self.assertEqual(metadata2.max_vel, 3.2)
        self.assertEqual(metadata2.min_vel, 2)

    def test_parse_static_solution_single_segment(self):
        # write sample single segment
        with open(self.solution_filepath, "w") as f:
            f.write(SINGLE_SEGMENT)

        # parse solution
        solution = Solution.from_static_solution_file(
            segments_data=SINGLE_SEGMENTS_DATA, solution_filepath=self.solution_filepath
        )

        # check values
        self.assertEqual(solution.num_segments, 1)
        self.assertEqual(solution.num_segments, len(solution.segments))
        segment = solution.segments[0]
        self.assertListEqual(list(segment.data.lat), [-32.851952, -32.691689])
        self.assertListEqual(list(segment.data.lon), [-72.503555, -72.481743])
        self.assertListEqual(list(segment.data.depth), [0.716800, 0.716800])
        self.assertListEqual(list(segment.data.slip), [126.725410, 144.565353])
        self.assertListEqual(list(segment.data.rake), [89.414207, 91.090263])
        self.assertListEqual(list(segment.data.strike), [6.613912, 6.613912])
        self.assertListEqual(list(segment.data.dip), [19.280828, 19.280828])
        self.assertListEqual(list(segment.data.trup), [84.339615, 75.787537])
        self.assertListEqual(list(segment.data.trise), [13.5, 9])
        self.assertListEqual(list(segment.data.tfall), [1.5, 1.5])
        self.assertListEqual(list(segment.data.moment), [0.119633e27, 0.136475e27])
        self.assertListEqual(
            segment.polygon,
            [
                (-72.5035553, -32.8519516, 0.716799974),
                (-71.2028580, -32.9804459, 44.0831985),
                (-70.7055817, -29.3264771, 44.0831985),
                (-72.0062790, -29.1979809, 0.716799974),
                (-72.5035553, -32.8519516, 0.716799974),
            ],
        )
        self.assertEqual(segment.nx, 23)
        self.assertEqual(segment.dx, 17.93)
        self.assertEqual(segment.ny, 9)
        self.assertEqual(segment.dy, 14.92)
        metadata = segment.metadata
        self.assertEqual(metadata.delay_segment, 0)
        self.assertEqual(metadata.delta_dip, 14)
        self.assertEqual(metadata.delta_strike, 17)
        self.assertEqual(metadata.dip, 19)
        self.assertEqual(metadata.dip_subfaults, 9)
        self.assertEqual(metadata.hyp_dip, 5)
        self.assertEqual(metadata.hyp_stk, 9)
        self.assertEqual(metadata.neighbours, [])
        self.assertEqual(metadata.rake, 109)
        self.assertEqual(metadata.rupture_vel, 2.5)
        self.assertEqual(metadata.stk_subfaults, 23)
        self.assertEqual(metadata.strike, 6)
        self.assertEqual(metadata.max_vel, 3.125)
        self.assertEqual(metadata.min_vel, 1)
