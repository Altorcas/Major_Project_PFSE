"""
A module for checkinG base plates
"""

from dataclasses import dataclass
from math import pi, sqrt
from handcalcs.decorator import handcalc


@dataclass
class baseplate:
    """
    A data type to represent the section of the column according to the Eurocode nomenclature
    and to represent base plate dimensions and characteristics

    h - height
    b - width
    tw - thickness of the web
    tf - thickness of the flange
    fyk - steel grade
    phi0 - factor of safety for steel
    d_edge - distance from the edge of the flange to the edge of the base plate
    H - dimension of the base plate paralel to the section height
    B - dimension of the base plate paralel to the section width
    t - thickness of the plate
    fck - characteristic compression resistance of the concrete
    phic - factor of safety for concrete
    h_stiff - stiffener height
    t_stiff - stiffener thickness
    fyk_stiff - stiffener steel grade
    D_bolt - bolt diameter
    L_bolt - bolt length
    d_bolt - distance from the base plate edge to the center of the bolt
    fyk_bolt - Yield strength for bolts
    phi_bolt - factor of safety for bolt steel
    E_bolt - Young modulus of bolts

    This is only valid for double T sections

    """
    h: float
    b: float
    tw: float
    tf: float
    fyk: float
    phi0: float
    d_edge: float
    H: float
    B: float
    t: float   
    fck: float
    phic: float
    h_stiff: float
    t_stiff: float
    fyk_stiff: float
    D_bolt: float
    L_bolt: float
    d_bolt: float
    fyk_bolt: float
    phi_bolt: float
    E_bolt: float

    
    def additional_bearing_width(self) -> float:
        """
        Returns the additional bearing with according to EC3 1-8 clause 6.2.5 (4)
        """

        c = self.t*sqrt((self.fyk/self.phi0)/(3*self.fck/self.phic))

        return c
    
    
    def effective_area(self) -> float:
        """
        Returns the effective area of the section with the base plate according to EC3 1-8 clause 6.2.5
        """
        c = self.additional_bearing_width()

        if self.H < self.h:
            return print("Section dimension h cannot be greater than base plate dimension H")
        
        if self.B < self.b:
            return print("Section dimension b cannot be greater than base plate dimension B")
            
        if c > (self.H/2-self.h/2) and c > (self.B/2-self.b/2):
            print("Adittional bearing width exceed H and B")
            return self.B*self.H-4*((self.B/2-c-self.tw/2)*(self.h/2-self.tf-c))
        if c > (self.H/2-self.h/2) and c < (self.B/2-self.b/2):
            print("Adittional bearing width exceed H but not B")
            return (self.tw+2*c)*(self.h-2*c-2*self.tf)+2*(self.b+2*c)*(c+self.tf+self.H/2-self.h/2)
        if c < (self.H/2-self.h/2) and c > (self.B/2-self.b/2):
            print("Adittional bearing width exceed B but not H")
            return (self.tw+2*c)*(self.h-2*c-2*self.tf)+2*(self.tf+2*c)*self.B
        if c < (self.H/2-self.h/2) and c < (self.B/2-self.b/2):
            print("Adittional bearing width dooes not exceed H nor B")
            return (2*c+self.tf)*(2*c+self.b)*2+(2*c+self.tw)*(self.h-2*self.tf-2*c)
        
    def Wel_effective_baseplate(self) -> float:
        """
        Calculates the effective elastic section modulus of the base plate, considering that the moment is applied around the
        strong axis of the steel section.
        """
        c = self.additional_bearing_width()

        if self.H < self.h:
            return print("Section dimension h cannot be greater than base plate dimension H")
        
        if self.B < self.b:
            return print("Section dimension b cannot be greater than base plate dimension B")

        if c > (self.H/2-self.h/2) and c > (self.B/2-self.b/2):

            second_moment_of_area = 1/12*(self.B*(self.H**3)-(self.B-(2*c+self.tw))*(self.H-2*(c+self.tf+self.H/2-self.h/2))**3)
            Wel = 2*second_moment_of_area/self.H

            return Wel

        if c > (self.H/2-self.h/2) and c < (self.B/2-self.b/2):

            second_moment_of_area = 1/12*((self.b+2*c)*self.H**3-((self.b+2*c)-(2*c+self.tw))*(self.H-2*(c+self.tf+self.H/2-self.h/2))**3)
            Wel = 2*second_moment_of_area/self.H

            return Wel
        
        if c < (self.H/2-self.h/2) and c > (self.B/2-self.b/2):

            second_moment_of_area = 1/12*((self.B*(self.h+2*c)**3)-((self.B-(2*c+self.tw)))*((self.h+2*c)-2*(2*c+self.tf))**3)
            Wel = 2*second_moment_of_area/(self.h+2*c)

            return Wel
        
        if c < (self.H/2-self.h/2) and c < (self.B/2-self.b/2):

            second_moment_of_area = 1/12*(((self.b+2*c)*(self.h+2*c)**3)-(((self.b+2*c)-(2*c+self.tw)))*((self.h+2*c)-2*(2*c+self.tf))**3)
            Wel = 2*second_moment_of_area/(self.h+2*c)

            return Wel
        

@handcalc()
def calc_additional_bearing_width(t, fyk, phi0, fck, phic) -> float:
    """
    Returns the effective area of the section with the base plate according to EC3 1-8 clause 6.2.5
    """
    c = (t*sqrt((fyk/phi0)/(3*fck/phic)))
       
    return c
            
        

def calc_effective_area(h, b, tw, tf, fyk, phi0, H, B, t, fck, phic) -> tuple[str, float]:
    """
    Returns the LaTeX representation of the effective area calculation formula
    and its numeric value according to EC3 1-8 clause 6.2.5.
    """
    c = t * sqrt((fyk / phi0) / (3 * fck / phic))
    formula_latex = ""
    A_eff = 0

    if c > (H / 2 - h / 2) and c > (B / 2 - b / 2):
        A_eff = B * H - 4 * ((B / 2 - c - tw / 2) * (h / 2 - tf - c))
        formula_latex = r"A_{eff} = B \cdot H - 4 \cdot \left( \left( \frac{B}{2} - c - \frac{tw}{2} \right) \cdot \left( \frac{h}{2} - tf - c \right) \right)"
    elif c > (H / 2 - h / 2) and c < (B / 2 - b / 2):
        A_eff = (tw + 2 * c) * (h - 2 * c - 2 * tf) + 2 * (b + 2 * c) * (c + tf + H / 2 - h / 2)
        formula_latex = r"A_{eff} = (tw + 2c) \cdot (h - 2c - 2tf) + 2 \cdot (b + 2c) \cdot (c + tf + \frac{H}{2} - \frac{h}{2})"
    elif c < (H / 2 - h / 2) and c > (B / 2 - b / 2):
        A_eff = (tw + 2 * c) * (h - 2 * c - 2 * tf) + 2 * (tf + 2 * c) * B
        formula_latex = r"A_{eff} = (tw + 2c) \cdot (h - 2c - 2tf) + 2 \cdot (tf + 2c) \cdot B"
    elif c < (H / 2 - h / 2) and c < (B / 2 - b / 2):
        A_eff = (2 * c + tf) * (2 * c + b) * 2 + (2 * c + tw) * (h - 2 * tf - 2 * c)
        formula_latex = r"A_{eff} = (2c + tf) \cdot (2c + b) \cdot 2 + (2c + tw) \cdot (h - 2tf - 2c)"

    return formula_latex, A_eff

@handcalc()
def ratios(M_Ed, N_Ed, A_eff, W_el) -> tuple:
    DesignValueRatio = M_Ed/N_Ed
    ResistanceRatio = W_el/A_eff
    return [DesignValueRatio, ResistanceRatio]

@handcalc()
def joint_design_bearing_strength(fck,phic) -> float:
    f_jd = fck/phic
    return f_jd

@handcalc()
def check_compression(NEd, Ap, MEd, Wp) -> float:
    Stress = NEd/Ap + MEd/Wp
    return Stress

@handcalc()
def bearing_width(t_stiff, b, B, c):
    b_bearing= 2*t_stiff+min(2*c,b)+min(2*c,B-b-2*t_stiff)
    return b_bearing

def compression_width(H, d_bolt, MEd, NEd, bearing_b, fck, phic):
    """
    Solves a quadratic ecuation to obtain block compression width
    """
    a = 1
    b = 2 * (-H+d_bolt)
    c = 2 * ((MEd + NEd*(H/2-d_bolt))/(bearing_b*fck/phic))

    discriminant = b**2-4*a*c

    if discriminant > 0:
        y1 = (-b + sqrt(discriminant)) / (2*a)
        y2 = (-b - sqrt(discriminant)) / (2*a)
        return y1, y2
    elif discriminant == 0:
        y = -b / (2*a)
        return y, y  # Both roots are the same
    else:
        real_part = -b / (2*a)
        imaginary_part = sqrt(abs(discriminant)) / (2*a)
        y1 = complex(real_part, imaginary_part)
        y2 = complex(real_part, -imaginary_part)
        return y1, y2
        
def bolt_stress(NEd, bearing_b, y, fck, phic) -> float:
    """
    Returns the base plate tenssion value
    """
    Td = bearing_b*y*fck/phic-NEd
    
    return Td

@handcalc()
def bp_stiff_gravity_center(h_stiff, t_stiff, t, B):
    """
    Returns the new stiffener + base plate gravity center
    """
    y_g = ((2*h_stiff*t_stiff)*(t+h_stiff/2)+B*t*t/2)/(2*h_stiff*t_stiff+B*t)

    return y_g

@handcalc(override = 'long')
def bp_stiff_momentinertia(h_stiff, t_stiff, t, B, y_g):
    """
    Returns the new stiffener + base plate moment of inertia
    """

    I_stiff = 1/12*t_stiff*h_stiff**3
    I_bp = 1/12*B*t**3
    d1 = 2*t_stiff*h_stiff*(t + h_stiff/2 - y_g)**2
    d2 = B*t*(y_g-t/2)**2

    I = 2*I_stiff + I_bp + d1 + d2

    return I

@handcalc(override = 'long')
def bp_stiff_wel(I, y_g, h_stiff, t):
    """
    Returns the new stiffener + base plate Well
    """

    W_el = min(I/y_g, I/(h_stiff+t-y_g))

    return W_el

def calculation_moment(fck, phic, bearing_b, y, d_edge, Td, d_bolt):
    """
    Returns the calculation bending moment to check base plate under compression
    """
    if y <= d_edge:
        MEd_calc_compression = fck/phic*y*(d_edge-y/2)

    if y > d_edge:
        MEd_calc_compression = fck/phic*bearing_b*(d_edge**2)/2

    MEd_calc_tenssion = Td*(d_edge - d_bolt)

    Med_calc = max(MEd_calc_compression, MEd_calc_tenssion)
        
    return Med_calc

@handcalc()
def uls_compression_bp_stiffener(M_Ed, W_el, f_yk, phi0):
    """
    Returns the ULS check for compression in the base plate (with stiffeners)
    """

    Ratio = M_Ed/(W_el*(f_yk/phi0))

    return Ratio

@handcalc()
def steel_stress_baseplate(y, H, d_edge, E_bolt) -> float:
    """
    Returns steel stress at the base plate to check if bolts are necessary or not.
    It is considered that unit stress for concrete is 0.0035.
    It is considered rectangular diagram for the equilibrium (y = 0.8x)
    """
    e_c = 0.0035
    x = y/0.8
    e_s = (H-d_edge-x)/x*e_c
    Steel_Stress = E_bolt*e_s

    return Steel_Stress

@handcalc()
def bolt_area(T_d, f_yk, phi_bolt) -> float:
    """
    Returns the required bolt surface per side
    """
    A_b = T_d/(f_yk/phi_bolt)

    return A_b

@handcalc()
def number_bolt_20(A_b):
    """
    Returns the number of bolts required (diameter 20mm)
    """
    n = A_b/(0.8*pi*20**2/4)

    return n

@handcalc()
def number_bolt_25(A_b):
    """
    Returns the number of bolts required (diameter 25mm)
    """
    n = A_b/(0.8*pi*25**2/4)
    
    return n

@handcalc()
def number_bolt_30(A_b):
    """
    Returns the number of bolts required (diameter 30mm)
    """
    n = A_b/(0.8*pi*30**2/4)
    
    return n

