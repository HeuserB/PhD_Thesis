import numpy as np
import itertools
import threading
import os
import logging
import lmfit

sin, cos, pi, ceil, sqrt, asin         =   np.sin, np.cos, np.pi, np.ceil, np.sqrt, np.arcsin 
logger = logging.getLogger(__name__)
epsilon = 1.0e-6  # for floating point comparison

import units

class Cell(object):
    """
    This is a cell object, able to calculate the volume and d-spacing according to formula from:

    http://geoweb3.princeton.edu/research/MineralPhy/xtalgeometry.pdf
    """
    lattices = ["cubic", "tetragonal", "hexagonal", "rhombohedral", "orthorhombic", "monoclinic", "triclinic"]
    types = {"P": "Primitive",
             "I": "Body centered",
             "F": "Face centered",
             "C": "Side centered",
             "R": "Rhombohedral"}

    def __init__(self, a=1, b=1, c=1, alpha=90, beta=90, gamma=90, lattice="triclinic", lattice_type="P"):
        """Constructor of the Cell class:

        Crystalographic units are Angstrom for distances and degrees for angles !

        :param a,b,c: unit cell length in Angstrom
        :param alpha, beta, gamma: unit cell angle in degrees
        :param lattice: "cubic", "tetragonal", "hexagonal", "rhombohedral", "orthorhombic", "monoclinic", "triclinic"
        :param lattice_type: P, I, F, C or R
        """
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.lattice = lattice if lattice in self.lattices else "triclinic"

        self._volume = None
        self.S11 = None
        self.S12 = None
        self.S13 = None
        self.S22 = None
        self.S23 = None
        self.selection_rules = []
        "contains a list of functions returning True(allowed)/False(forbiden)/None(unknown)"
        self._type = "P"
        self.set_type(lattice_type)

    def __repr__(self, *args, **kwargs):
        return "%s %s cell a=%.4f b=%.4f c=%.4f alpha=%.3f beta=%.3f gamma=%.3f" % \
            (self.types[self.type], self.lattice, self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    @classmethod
    def cubic(cls, a, lattice_type="P"):
        """Factory for cubic lattices

        :param a: unit cell length
        """
        a = float(a)
        self = cls(a, a, a, 90, 90, 90,
                   lattice="cubic", lattice_type=lattice_type)
        return self

    @classmethod
    def tetragonal(cls, a, c, lattice_type="P"):
        """Factory for tetragonal lattices

        :param a: unit cell length
        :param c: unit cell length
        """
        a = float(a)
        self = cls(a, a, float(c), 90, 90, 90,
                   lattice="tetragonal", lattice_type=lattice_type)
        return self

    @classmethod
    def orthorhombic(cls, a, b, c, lattice_type="P"):
        """Factory for orthorhombic lattices

        :param a: unit cell length
        :param b: unit cell length
        :param c: unit cell length
        """
        self = cls(float(a), float(b), float(c), 90, 90, 90,
                   lattice="orthorhombic", lattice_type=lattice_type)
        return self

    @classmethod
    def hexagonal(cls, a, c, lattice_type="P"):
        """Factory for hexagonal lattices

        :param a: unit cell length
        :param c: unit cell length
        """
        a = float(a)
        self = cls(a, a, float(c), 90, 90, 120,
                   lattice="hexagonal", lattice_type=lattice_type)
        return self

    @classmethod
    def monoclinic(cls, a, b, c, beta, lattice_type="P"):
        """Factory for hexagonal lattices

        :param a: unit cell length
        :param b: unit cell length
        :param c: unit cell length
        :param beta: unit cell angle
        """
        self = cls(float(a), float(b), float(c), 90, float(beta), 90,
                   lattice_type=lattice_type, lattice="monoclinic")
        return self

    @classmethod
    def rhombohedral(cls, a, alpha, lattice_type="P"):
        """Factory for hexagonal lattices

        :param a: unit cell length
        :param alpha: unit cell angle
        """
        a = float(a)
        alpha = float(a)
        self = cls(a, a, a, alpha, alpha, alpha,
                   lattice="rhombohedral", lattice_type=lattice_type)
        return self

    @classmethod
    def diamond(cls, a):
        """Factory for Diamond type FCC like Si and Ge

        :param a: unit cell length
        """
        self = cls.cubic(a, lattice_type="F")
        self.selection_rules.append(lambda h, k, l: not((h % 2 == 0) and (k % 2 == 0) and (l % 2 == 0) and ((h + k + l) % 4 != 0)))
        return self

    @property
    def volume(self):
        if self._volume is None:
            self._volume = self.a * self.b * self.c
            if self.lattice not in ["cubic", "tetragonal", "orthorhombic"]:
                cosa = cos(self.alpha * pi / 180.)
                cosb = cos(self.beta * pi / 180.)
                cosg = cos(self.gamma * pi / 180.)
                self._volume *= sqrt(1 - cosa ** 2 - cosb ** 2 - cosg ** 2 + 2 * cosa * cosb * cosg)
        return self._volume

    def get_type(self):
        return self._type

    def set_type(self, lattice_type):
        self._type = lattice_type if lattice_type in self.types else "P"
        self.selection_rules = [lambda h, k, l: not(h == 0 and k == 0 and l == 0)]
        if self._type == "I":
            self.selection_rules.append(lambda h, k, l: (h + k + l) % 2 == 0)
        if self._type == "F":
            self.selection_rules.append(lambda h, k, l: (h % 2 + k % 2 + l % 2) in (0, 3))
        if self._type == "R":
            self.selection_rules.append(lambda h, k, l: ((h - k + l) % 3 == 0))

    type = property(get_type, set_type)

    def d(self, hkl):
        """
        Calculate the actual d-spacing for a 3-tuple of integer representing a
        family of Miller plans

        :param hkl: 3-tuple of integers
        :return: the inter-planar distance
        """
        h, k, l = hkl
        if self.lattice in ["cubic", "tetragonal", "orthorhombic"]:
            invd2 = (h / self.a) ** 2 + (k / self.b) ** 2 + (l / self.c) ** 2
        else:
            if self.S11 is None:
                alpha = self.alpha * pi / 180.
                cosa = cos(alpha)
                sina = sin(alpha)
                beta = self.beta * pi / 180.
                cosb = cos(beta)
                sinb = sin(beta)
                gamma = self.gamma * pi / 180.
                cosg = cos(gamma)
                sing = sin(gamma)

                self.S11 = (self.b * self.c * sina) ** 2
                self.S22 = (self.a * self.c * sinb) ** 2
                self.S33 = (self.a * self.b * sing) ** 2
                self.S12 = self.a * self.b * self.c * self.c * (cosa * cosb - cosg)
                self.S23 = self.a * self.a * self.b * self.c * (cosb * cosg - cosa)
                self.S13 = self.a * self.b * self.b * self.c * (cosg * cosa - cosb)

            invd2 = (self.S11 * h * h +
                     self.S22 * k * k +
                     self.S33 * l * l +
                     2 * self.S12 * h * k +
                     2 * self.S23 * k * l +
                     2 * self.S13 * h * l)
            invd2 /= (self.volume) ** 2
        return sqrt(1 / invd2)

    def d_spacing(self, dmin=1.0):
        """Calculate all d-spacing down to dmin

        applies selection rules

        :param dmin: minimum value of spacing requested
        :return: dict d-spacing as string, list of tuple with Miller indices
                preceded with the numerical value
        """
        hmax = int(ceil(self.a / dmin))
        kmax = int(ceil(self.b / dmin))
        lmax = int(ceil(self.c / dmin))
        res = {}
        for hkl in itertools.product(range(-hmax, hmax + 1),
                                     range(-kmax, kmax + 1),
                                     range(-lmax, lmax + 1)):
            # Apply selection rule
            valid = True
            for rule in self.selection_rules:
                valid = rule(*hkl)
                if not valid:
                    break
            if not valid:
                continue

            d = self.d(hkl)
            strd = "%.8e" % d
            if d < dmin:
                continue
            if strd in res:
                res[strd].append(hkl)
            else:
                res[strd] = [d, hkl]
        return res

class Calibrant(object):
    """
    A calibrant is a reference compound where the d-spacing (interplanar distances)
    are known. They are expressed in Angstrom (in the file)
    """

    def __init__(self, filename=None, dSpacing=None, wavelength=None):
        object.__init__(self)
        self._filename = filename
        self._wavelength = wavelength
        self._sem = threading.Semaphore()
        self._2th = []
        if filename is not None:
            self._dSpacing = None
        elif dSpacing is None:
            self._dSpacing = []
        else:
            self._dSpacing = list(dSpacing)
        self._out_dSpacing = []
        if self._dSpacing and self._wavelength:
            self._calc_2th()

    def __eq__(self, other):
        """
        Test the equality with another object

        It only takes into acount the wavelength and dSpacing, not the
        filename.

        :param object other: Another object
        :rtype: bool
        """
        if other is None:
            return False
        if not isinstance(other, Calibrant):
            return False
        if self._wavelength != other._wavelength:
            return False
        if self.dSpacing != other.dSpacing:
            return False
        return True

    def __ne__(self, other):
        """
        Test the non-equality with another object

        It only takes into acount the wavelength and dSpacing, not the
        filename.

        :param object other: Another object
        :rtype: bool
        """
        return not (self == other)

    def __hash__(self):
        """
        Returns the hash of the object.

        It only takes into acount the wavelength and dSpacing, not the
        filename.

        :rtype: int
        """
        h = hash(self._wavelength)
        for d in self.dSpacing:
                h = h ^ hash(d)
        return h

    def __copy__(self):
        """
        Copy a calibrant

        :rtype: Calibrant
        """
        self._initialize()
        return Calibrant(filename=self._filename,
                         dSpacing=self._dSpacing + self._out_dSpacing,
                         wavelength=self._wavelength)

    def __repr__(self):
        name = "undefined"
        if self._filename:
            name = os.path.splitext(os.path.basename(self._filename))[0]
        name += " Calibrant "
        if len(self.dSpacing):
            name += "with %i reflections " % len(self._dSpacing)
        if self._wavelength:
            name += "at wavelength %s" % self._wavelength
        return name

    def get_filename(self):
        return self._filename

    filename = property(get_filename)

    def load_file(self, filename=None):
        with self._sem:
            self._load_file(filename)

    def _load_file(self, filename=None):
        if filename:
            self._filename = filename
        if not os.path.isfile(self._filename):
            logger.error("No such calibrant file: %s", self._filename)
            return
        self._filename = os.path.abspath(self._filename)
        self._dSpacing = np.unique(np.loadtxt(self._filename))
        self._dSpacing = list(self._dSpacing[-1::-1])  # reverse order
        # self._dSpacing.sort(reverse=True)
        if self._wavelength:
            self._calc_2th()

    def _initialize(self):
        """Initialize the object if expected."""
        if self._dSpacing is None:
            if self._filename:
                self._load_file()
            else:
                self._dSpacing = []

    def count_registered_dSpacing(self):
        """Count of registered dSpacing positons."""
        self._initialize()
        return len(self._dSpacing) + len(self._out_dSpacing)

    def save_dSpacing(self, filename=None):
        """
        save the d-spacing to a file

        """
        self._initialize()
        if (filename is None) and (self._filename is not None):
            filename = self._filename
        else:
            return
        with open(filename) as f:
            f.write("# %s Calibrant" % filename)
            for i in self.dSpacing:
                f.write("%s\n" % i)

    def get_dSpacing(self):
        self._initialize()
        return self._dSpacing

    def set_dSpacing(self, lst):
        self._dSpacing = list(lst)
        self._out_dSpacing = []
        self._filename = "Modified"
        if self._wavelength:
            self._calc_2th()

    dSpacing = property(get_dSpacing, set_dSpacing)

    def append_dSpacing(self, value):
        self._initialize()
        with self._sem:
            delta = [abs(value - v) / v for v in self._dSpacing if v is not None]
            if not delta or min(delta) > epsilon:
                self._dSpacing.append(value)
                self._dSpacing.sort(reverse=True)
                self._calc_2th()

    def append_2th(self, value):
        with self._sem:
            self._initialize()
            if value not in self._2th:
                self._2th.append(value)
                self._2th.sort()
                self._calc_dSpacing()

    def setWavelength_change2th(self, value=None):
        with self._sem:
            if value:
                self._wavelength = float(value)
                if self._wavelength < 1e-15 or self._wavelength > 1e-6:
                    logger.warning("This is an unlikely wavelength (in meter): %s", self._wavelength)
                self._calc_2th()

    def setWavelength_changeDs(self, value=None):
        """
        This is probably not a good idea, but who knows !
        """
        with self._sem:
            if value:
                self._wavelength = float(value)
                if self._wavelength < 1e-15 or self._wavelength > 1e-6:
                    logger.warning("This is an unlikely wavelength (in meter): %s", self._wavelength)
                self._calc_dSpacing()

    def set_wavelength(self, value=None):
        updated = False
        with self._sem:
            if self._wavelength is None:
                if value:
                    self._wavelength = float(value)
                    if (self._wavelength < 1e-15) or (self._wavelength > 1e-6):
                        logger.warning("This is an unlikely wavelength (in meter): %s", self._wavelength)
                    updated = True
            elif abs(self._wavelength - value) / self._wavelength > epsilon:
                logger.warning("Forbidden to change the wavelength once it is fixed !!!!")
                logger.warning("%s != %s, delta= %s", self._wavelength, value, self._wavelength - value)
        if updated:
            self._calc_2th()

    def get_wavelength(self):
        return self._wavelength

    wavelength = property(get_wavelength, set_wavelength)

    def _calc_2th(self):
        """Calculate the 2theta positions for all peaks"""
        self._initialize()
        if self._wavelength is None:
            logger.error("Cannot calculate 2theta angle without knowing wavelength")
            return
        tths = []
        dSpacing = self._dSpacing[:] + self._out_dSpacing  # explicit copy
        try:
            for ds in dSpacing:
                tth = 2.0 * asin(5.0e9 * self._wavelength / ds)
                tths.append(tth)
        except ValueError:
            size = len(tths)
            # remove dSpacing outside of 0..180
            self._dSpacing = dSpacing[:size]
            self._out_dSpacing = dSpacing[size:]
        else:
            self._dSpacing = dSpacing
            self._out_dSpacing = []
        self._2th = tths

    def _calc_dSpacing(self):
        if self._wavelength is None:
            logger.error("Cannot calculate 2theta angle without knowing wavelength")
            return
        self._dSpacing = [5.0e9 * self._wavelength / sin(tth / 2.0) for tth in self._2th]

    def get_2th(self):
        """Returns the 2theta positions for all peaks (cached)"""
        if not self._2th:
            self._initialize()
            if not self._dSpacing:
                logger.error("Not d-spacing for calibrant: %s", self)
            with self._sem:
                if not self._2th:
                    self._calc_2th()
        return self._2th

    def get_2th_index(self, angle, delta=None):
        """Returns the index in the 2theta angle index

        :param angle: expected angle in radians
        :param delta: precision on angle
        :return: 0-based index or None
        """
        if angle and angle in self._2th:
            return self._2th.index(angle)
        if delta:
            d2th = abs(np.array(self._2th) - angle)
            if d2th.min() < delta:
                return d2th.argmin()

    def get_max_wavelength(self, index=None):
        """Calculate the maximum wavelength assuming the ring at index is visible

        Bragg's law says: $\\lambda = 2d sin(\\theta)$
        So at 180° $\\lambda = 2d$

        :param index: Ring number, otherwise assumes all rings are visible
        :return: the maximum visible wavelength
        """
        dSpacing = self._dSpacing[:] + self._out_dSpacing  # get all rings
        if index is None:
            index = len(dSpacing) - 1
        if index >= len(dSpacing):
            raise IndexError("There are not than many (%s) rings indices in this calibrant" % (index))
        return dSpacing[index] * 2e-10

    def get_peaks(self, unit="2th_deg"):
        """Calculate the peak position as
        :return: numpy array (unlike other methods which return lists)
        """
        unit = units.to_unit(unit)
        scale = unit.scale
        name = unit.name
        size = len(self.get_2th())
        if name.startswith("2th"):
            values = np.array(self.get_2th())
        elif name.startswith("q"):
            values = 20.0 * pi / np.array(self.get_dSpacing()[:size])
        else:
            raise ValueError("Only 2\theta and *q* units are supported for now")

        return values * scale

    def fake_calibration_image(self, ai, shape=None, Imax=1.0,
                               U=0, V=0, W=0.0001):
        """
        Generates a fake calibration image from an azimuthal integrator

        :param ai: azimuthal integrator
        :param Imax: maximum intensity of rings
        :param U, V, W: width of the peak from Caglioti's law (FWHM^2 = Utan(th)^2 + Vtan(th) + W)
        """
        if shape is None:
            if ai.detector.shape:
                shape = ai.detector.shape
            elif ai.detector.max_shape:
                shape = ai.detector.max_shape
        if shape is None:
            raise RuntimeError("No shape available")
        if (self.wavelength is None) and (ai._wavelength is not None):
            self.wavelength = ai.wavelength
        elif (self.wavelength is None) and (ai._wavelength is None):
            raise RuntimeError("Wavelength needed to calculate 2theta position")
        elif (self.wavelength is not None) and (ai._wavelength is not None) and\
                abs(self.wavelength - ai.wavelength) > 1e-15:
            logger.warning("Mismatch between wavelength for calibrant (%s) and azimutal integrator (%s)",
                           self.wavelength, ai.wavelength)
        tth = ai.twoThetaArray(shape)
        tth_min = tth.min()
        tth_max = tth.max()
        dim = int(np.sqrt(shape[0] * shape[0] + shape[1] * shape[1]))
        tth_1d = np.linspace(tth_min, tth_max, dim)
        tanth = np.tan(tth_1d / 2.0)
        fwhm2 = U * tanth ** 2 + V * tanth + W
        sigma2 = fwhm2 / (8.0 * np.log(2.0))
        signal = np.zeros_like(sigma2)
        for t in self.get_2th():
            if t >= tth_max:
                break
            else:
                signal += Imax * np.exp(-(tth_1d - t) ** 2 / (2.0 * sigma2))
        res = ai.calcfrom1d(tth_1d, signal, shape=shape, mask=ai.mask,
                            dim1_unit='2th_rad', correctSolidAngle=True)
        return res

    def __getnewargs_ex__(self):
        return (self._filename, self._dSpacing, self._wavelength), {}

    def __getstate__(self):
        state_blacklist = ('_sem',)
        state = self.__dict__.copy()
        for key in state_blacklist:
            if key in state:
                del state[key]
        return state

    def __setstate__(self, state):
        for statekey, statevalue in state.items():
            setattr(self, statekey, statevalue)
        self._sem = threading.Semaphore()

class ImageGenerator(object):
    """
    An ImageGenerator is a class which takes a calibration, detector and calbibrant
    and simulates 2D images from this
    """
    def __init__(self, poni_file, detector_dict, calibrant_file, mask, wavelength,scale=1e4):
        self.wavelength =   wavelength
        self.mask   =   mask
        self.scale  =   scale
        self.calibrant  =   Calibrant()
        self.calibrant.load_file(calibrant_file)
        self.calibrant.set_wavelength(wavelength)
        self.ai =   pyFAI.load(poni_file_SACLA)
        self.detector   =   pyFAI.detectors.Detector(pixel1=detector_dict["pixel1"], pixel2=detector_dict["pixel2"], max_shape=detector_dict['max_shape'])
        self.ai.detector =   self.detector
        self.Image2d    =   None
        self.Noise2d    =   None
        self.__setImage2d__()
 
    def simulate2d(self):
        fake_image          =   self.calibrant.fake_calibration_image(self.ai)
        fake_image          =   np.flipud(fake_image)
        return fake_image
    
    def computeLineout(self):
        if self.Image2d is None:
            self.__setImage2d__()
        lineout =   self.ai.integrate1d(np.flipud(self.Image2d),npt=1000,method='csr',unit='2th_deg',correctSolidAngle=True,\
                        polarization_factor=0.99,mask=self.mask,error_model = "azimuthal")
        return np.array(lineout)
    
    def computeNoise(self,noisetype):
        """ Calculate noise on top of the image assuming a certain distribution of the statistics
        noisetype = 'Poisson', variance: $\lambda = \sqrt{sigma}$
        """
        if noisetype=="Poisson":
            sigma    =   np.sqrt(self.Image2d)
            #self.Noise2d   =   np.random.normal(np.zeros_like(self.Image2d), sigma)
            self.Noise2d   =   np.random.poisson(self.Image2d)
        else:
            logger.error(f"Unknown noise type {noisetype}!")
            return
    
    def noisyLineout(self, noisetype):
        if self.Noise2d is None:
             self.computeNoise(noisetype)
        noisyLineout    =   self.ai.integrate1d(np.flipud(self.Image2d + self.Noise2d),npt=1000,method='csr',unit='2th_deg',correctSolidAngle=True,\
                        polarization_factor=0.99,mask=self.mask,error_model = "azimuthal")
        return np.array(noisyLineout)
    
    def makeNoiseStatistics(self, noisetype, n):
        lineout =   self.computeLineout()
        result  =   np.zeros([n,lineout.shape[1]])
        errors  =   np.zeros([n,lineout.shape[1]])
        if noisetype=="Poisson":
            for id in range(n):
                self.computeNoise("Poisson")
                tmp =   self.noisyLineout("Poisson")
                result[id]   =   tmp[1]
                errors[id]   =   tmp[2]
            return lineout, result, errors
                
    def __setImage2d__(self):
        self.Image2d    =   np.array(self.simulate2d() * self.scale,dtype=int)
        self.Noise2d    =   np.zeros_like(self.Image2d)

def fit_gaussian(data):
    x, y    =   data[0], data[1]
    gmodel  =   lmfit.Model(lmfit.models.gaussian)
    params  =   gmodel.make_params()
    params  =   gmodel.make_params(cen=0.0, amp=1, wid=0.25)
    result  =   gmodel.fit(y, params, x=x)
    return result

def fit_lorentzian(data):
    x, y    =   data[0], data[1]
    gmodel  =   lmfit.Model(lmfit.models.lorentzian)
    params  =   gmodel.make_params()
    params  =   gmodel.make_params(cen=0.0, amp=1, wid=0.25)
    result  =   gmodel.fit(y, params, x=x)
    return result


if __name__=="__main__":
    import pyFAI
    import imageio

    detector_dict   =   {'pixel1':50E-6, 'pixel2':50E-6, 'max_shape':[1548,2064]}
    calibrant_file  =   '/home/benjamin/.virtualenvs/phd_thesis/lib/python3.8/site-packages/pyFAI/resources/calibration/CeO2.D'
    poni_file_SACLA = "/home/benjamin//Nextcloud/Work/W_PhD/W_PhD_Git/PhD_Git/PhD_Thesis/.data_SACLA/poni/Q_SACLA.poni"
    mask_SACLA  = "../../.data_SACLA/mask/SACLA_extensive.mask"
    SACLA_mask  = np.array(imageio.imread(mask_SACLA),dtype=int);
    wavelength  =   1.2402999999999998e-10
    SACLA_generator =   ImageGenerator(poni_file_SACLA,detector_dict,calibrant_file,SACLA_mask,wavelength)
    
    lineout, result, errors =   SACLA_generator.makeNoiseStatistics('Poisson',100)
    np.savetxt('../../.data_SACLA/simulations/lineout.txt',lineout)
    np.savetxt('../../.data_SACLA/simulations/noisy_lineouts.txt',result)
    np.savetxt('../../.data_SACLA/simulations/error_pyFAI.txt',errors)

    lineout =   np.loadtxt('../../.data_SACLA/simulations/lineout.txt')
    result  =   np.loadtxt('../../.data_SACLA/simulations/noisy_lineouts.txt')
    errors  =   np.loadtxt('../../.data_SACLA/simulations/error_pyFAI.txt')

    import matplotlib.pyplot as plt
    fig, ax =   plt.subplots(1,3)
    residuals = result - lineout[1]
    checkpoints =   np.array([400,720])

    ax[0].plot(lineout[0],lineout[1],c='k')
    for i in range(2):
        ax[0].plot(lineout[0],result[i])
    
    print(lineout[0][720])
    for pointID, point in enumerate(checkpoints):
        hist, edges    =   np.histogram(residuals[:,point],bins=50)
        edges   =   (edges[:-1] + edges[1:])/2.
        fit_result_gaussian =   fit_gaussian([edges,hist])
        #fit_result_lorentzian =   fit_lorentzian([edges,hist])
        ax[pointID+1].hist(residuals[:,point],bins=50)
        ax[pointID+1].scatter(edges, hist,c='k')
        ax[pointID+1].plot(edges, fit_result_gaussian.best_fit, '--', label='Gaussian fit',c='r')
        #ax[pointID+1].plot(edges, fit_result_lorentzian.best_fit, '--', label='Lorentzian fit',c='g')
        ax[pointID+1].legend()
        print(fit_result_gaussian.fit_report())
        ax[pointID+1].text(0.1,0.9,r'Gaussian: $\sigma=%.3f$\npiFAI error: $\partial_I={%.3f}$' %(fit_result_gaussian.values['sigma'],lineout[2,point]),transform=ax[pointID+1].transAxes)
    #ax[2].hist(residuals[:,720],bins=50)
    #ax[0].imshow(noisyImage)
    #ax[1].imshow(noise_image)
    #ax[2].plot(lineout_noise[0],lineout_noise[2],markersize=0.1, label="Error estimate from pyFAI")
    #ax[2].plot(lineout[0],sqrt(lineout[1]),label=r"$\sqrt{I}$")
    #ax[2].legend()
    plt.show()
    print(np.mean(SACLA_generator.Noise2d))