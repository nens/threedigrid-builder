class Lines:
    id = None # Array
    kcu = None # Array
    lik = None # Array  2D only for now.
    lim = None # Array 2D only for now.
    lin = None # Array 2D only for now.
    line = None # Array
    content_pk = None # Array
    content_type = None # Array
    line_geometries = None # Array
    dpumax = None # Array
    flod = None # Array Obstacle height for 2D for now.
    flou = None # Array Obstacle height for 2D for now.
    cross1 = None # Array cross-section definition id. # TODO: Discuss wether this is more table data or grid admin data.
    cross2 = None # Array cross-section definition id. TODO: Discuss wether this is more table data or grid admin data.
    cross_weight = None # Array interpolation weight for cross1. TODO: Discuss wether this is more table data or grid admin data.
    ds1d = None # Array Arc length
    fi1d = None # Array Angle of channel fork ## TODO: Get details of angle calc.
    ### TODO: Extend with extra info from sqlite such as zoom category etc.