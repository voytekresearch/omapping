"""Utility and helper functions for the testing of OM."""

import os
import pkg_resources as pkg

from om.meg.single import MegData
from om.meg.group import GroupMegData
from om.maps.anat import MapCompAnat
from om.core.db import OMDB
from om.core.osc import Osc

#############################################################################################
#############################################################################################
#############################################################################################

class TestDB(OMDB):
    """   """

    def __init__(self):

        # Initialize from OMDB object
        OMDB.__init__(self)

        # Set up the base path to tests data
        test_dat_path = pkg.resource_filename(__name__, 'data')
        # Note: the following also works:
        #   test_dat_path = pkg.resource_filename('om', 'tests/data')

        # Set base paths, and add 'Other' path, for test specific files
        self.internal_path = os.path.join(test_dat_path, 'Internal')
        self.external_path = os.path.join(test_dat_path, 'External')
        self.other_path = os.path.join(test_dat_path, 'Other')

        # Generate test paths
        self.gen_paths()

        # Add relevant paths to 'Other'
        self.csvs_path = os.path.join(self.other_path, 'csvs')

#############################################################################################
#############################################################################################
#############################################################################################

def load_test_meg_subj(sub, all_oscs=False, bands_vertex=False):
    """Loads a test subject of MD_SING data."""

    tdb = TestDB()
    osc = Osc(default=True)

    dat = MegData(tdb, '', osc)

    dat.import_foof(sub, get_demo=False, load_type='pickle')

    if all_oscs:
        dat.all_oscs()

    if bands_vertex:
        dat.osc_bands_vertex()

    return dat

def load_test_meg_gr(bands_vertex=False, all_osc=False, peaks=False, calc_maps=False, vertex_osc=False):
    """Loads a test group object of MD_GR data."""

    tdb = TestDB()
    osc = Osc(default=True)

    meg_group = GroupMegData(tdb, osc)

    subjs = ['test_v5', 'test_v5']

    for s in subjs:

        meg_subj = load_test_meg_subj(s)

        if bands_vertex:
            meg_subj.osc_bands_vertex()

        if all_osc:
            meg_subj.all_oscs()

        if peaks:
            meg_subj.peak_freq(dat='all')

        meg_group.add_subject(meg_subj, add_all_oscs=all_osc,
                              add_vertex_bands=bands_vertex,
                              add_peak_freqs=peaks,
                              add_vertex_oscs=vertex_osc)

    if calc_maps:
        meg_group.osc_prob()
        meg_group.osc_score()

    return meg_group

def load_test_anat(load_meg=False, load_scout=False, load_anat=False,
                   align=False, convert=False, calc_meg=False):
    """   """

    tdb = TestDB()

    map_comp = MapCompAnat(tdb)

    if load_meg:
        map_comp.load_meg_maps('test_meg')
    if load_scout:
        map_comp.load_elec_rois('test_scout.mat')
    if load_anat:
        map_comp.load_anat_maps('test_anat.mat', 'test')

    if align:
        map_comp.align_rois()
    if convert:
        map_comp.conv_meg_rois()
    if calc_meg:
        map_comp.calc_meg_con()

    return map_comp
