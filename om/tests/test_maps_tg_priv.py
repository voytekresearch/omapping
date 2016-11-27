

####################################################################################
################ TESTS - OMEGAMAPPIN - CL_MC_TG - PRIVATE FUNCTIONS ################
####################################################################################

def test_get_map_names_terms():
    """   """

    db = OMDB()

    names_file = '00-ns_terms.csv'
    names = mc._get_map_names(names_file, os.path.join(db.maps_path, 'Terms'))

    assert names
    assert type(names) == ListType
    assert type(names[0]) == StringType

def test_get_map_names_genes():
    """   """

    db = OMDB()

    names_file = '00-real_gene_names.csv'
    names = mc._get_map_names(names_file, os.path.join(db.maps_path, 'Genes'))

    assert names
    assert type(names) is ListType
    assert type(names[0]) is StringType

def test_get_gene_files():
    """   """

    db = OMDB()

    genes_files_path = mc._get_gene_files('sub1', db)

    assert genes_files_path
    assert type(genes_files_path) is ListType
    for path in genes_files_path:
        assert os.path.exists(path)

def test_avg_csv_files():
    """   """

    tdb = TDB()

    f_in = [os.path.join(tdb.csvs_path, 'test1.csv'),
            os.path.join(tdb.csvs_path, 'test2.csv')]

    f_out = os.path.join(tdb.csvs_path, 'test_out.csv')

    mc._avg_csv_files(f_in, f_out)

    assert os.path.exists(f_out)

    exp = [[1.5, 1.5, 1.5, 1.5], [2.5, 2.5, 2.5, 2.5]]

    f = open(f_out)
    reader = csv.reader(f)

    for ind, row in enumerate(reader):
        row = [float(i) for i in row]
        assert row == exp[ind]

    f.close()
    os.remove(f_out)

    mc._avg_csv_files(f_in, f_out, 'median')
    assert os.path.exists(f_out)
    os.remove(f_out)

def test_init_stat_dict():
    """   """

    bands = ['a', 'b', 'c']
    d = mc._init_stat_dict(bands)

    assert d
    bands.append('Slopes')
    assert set(d.keys()) == set(bands)

def test_make_list():
    """   """

    df = pd.DataFrame(np.array([[1, 2], [1, 2]]))

    l = mc._make_list(df)

    assert len(l) == 2
    assert np.array_equal(l[0], np.array([1, 1]))
    assert np.array_equal(l[1], np.array([2, 2]))

def test_pull_out_results():
    """   """

    dat = [(1, 2), (1, 2)]

    out_1, out_2 = mc._pull_out_results(dat)

    assert np.array_equal(out_1, np.array([1, 1]))
    assert np.array_equal(out_2, np.array([2, 2]))
