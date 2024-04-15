import base_plates as bp

def test_calc_additional_bearing_width():
    c = 50.35
    t = 22
    fyk = 275
    phi0 = 1.05
    fck = 25
    phic = 1.5
    assert round(bp.calc_additional_bearing_width(t, fyk, phi0, fck, phic)[1],2) == c

def test_calc_effective_area():
    Aeff = 93876.46
    h = 400
    b = 180
    tw = 8.6
    tf = 13.5
    fyk = 275
    phi0 = 1.05
    H = 600
    B = 400
    t = 22
    fck = 25
    phic = 1.5

    assert round(bp.calc_effective_area(h, b, tw, tf, fyk, phi0, H, B, t, fck, phic)[1],2) == Aeff   


def test_ratios():
    Ratio1 = 2614.61
    Ratio2 = 112.67
    Med = 198261000
    Ned = 75828
    Wel = 10576915.65
    Aeff = 93876.46

    assert round(bp.ratios(Med, Ned, Aeff, Wel)[1][0],2) == Ratio1
    assert round(bp.ratios(Med, Ned, Aeff, Wel)[1][1],2) == Ratio2

def test_joint_design_bearing_strength():
    fjd = 16.67
    fck = 25
    phic = 1.5

    assert round(bp.joint_design_bearing_strength(fck, phic)[1],2) == fjd
    
def test_check_compression():
    Stress = 19.55
    Med = 198261000
    Ned = 75828
    Wel = 10576915.65
    Aeff = 93876.46

    assert round(bp.check_compression(Ned, Aeff, Med, Wel)[1],2) == Stress    

def test_bearing_width():
    b_bearing = 221.40
    t_stiff = 10
    b = 180
    B = 400
    c = 50.35

    assert round(bp.bearing_width(t_stiff, b, B, c)[1],2) == b_bearing 

def test_compression_width():
    y1 = 979.84
    y2 = 120.16
    H = 600
    d_bolt = 50
    Med = 198261000
    Ned = 75828
    b_bearing = 221.40
    fck = 25
    phic = 1.5

    assert round(bp.compression_width(H, d_bolt, Med, Ned, b_bearing, fck, phic)[0],2) == y1
    assert round(bp.compression_width(H, d_bolt, Med, Ned, b_bearing, fck, phic)[1],2) == y2

def test_bolot_stress():
    Td = 367533.72
    Ned = 75828
    b_bearing = 221.40
    y = 120.15222661026536
    fck = 25
    phic = 1.5

    assert round(bp.bolt_stress(Ned, b_bearing, y, fck, phic),2) == Td

def test_bp_stiff_gravity_center():
    yg = 32.86
    h_stiff = 150
    t_stiff = 10
    t = 22
    B = 400

    assert round(bp.bp_stiff_gravity_center(h_stiff, t_stiff, t, B)[1],2) == yg

def test_bp_stiff_momentinertia():
    Ig = 22526916.6
    h_stiff = 150
    t_stiff = 10
    t = 22
    B = 400
    yg = 32.86

    assert round(bp.bp_stiff_momentinertia(h_stiff, t_stiff, t, B, yg)[1],1) == Ig

def test_bp_stiff_wel():
    Wel = 161901
    I = 22526916
    yg = 32.86
    h_stiff = 150
    t = 22

    assert round(bp.bp_stiff_wel(I, yg, h_stiff, t)[1],0) == Wel

def test_calculation_moment():
    Med_calc = 18450000
    y = 120.15222661026536
    fck = 25
    phic = 1.5
    b_bearing = 221.40
    d_edge = 100
    Td = 367533.72
    d_bolt = 50

    assert round(bp.calculation_moment(fck ,phic, b_bearing, y, d_edge, Td, d_bolt),1) == Med_calc

def test_uls_compression_bp_stiffener():
    ratio = 0.435
    Med_calc = 18450000
    Wel = 161901
    fyk = 275
    phi0 = 1.05

    assert round(bp.uls_compression_bp_stiffener(Med_calc, Wel, fyk, phi0)[1],3) == ratio

def test_steel_stress_baseplate():
    stress = 1956.64
    y = 120.15
    H = 600
    d_edge = 50
    E_bolt = 210000

    assert round(bp.steel_stress_baseplate(y, H, d_edge, E_bolt)[1],2) == stress

def test_bolt_area():
    A = 1536.96
    Td = 367533.72
    fyk = 275
    phi_bolt = 1.15

    assert round(bp.bolt_area(Td, fyk, phi_bolt)[1],2) == A    

def test_number_bolt_20():
    A = 1536.96
    n = 6.12

    assert round(bp.number_bolt_20(A)[1],2) == 6.12 

def test_number_bolt_25():
    A = 1536.96
    n = 6.12

    assert round(bp.number_bolt_25(A)[1],2) == 3.91 

def test_number_bolt_30():
    A = 1536.96
    n = 6.12

    assert round(bp.number_bolt_30(A)[1],2) == 2.72 

    

