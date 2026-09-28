# Magnetosphere / ionosphere current explorer

Published independently of BELLA at https://juha.no/currents/.
BELLA's dark canvas, orbit interaction and compact panel are the visual reference;
its measurement projection is not used. Three.js 0.186.1 is vendored with license.

Regenerate data with `conda run -n base python build_traces.py` (numpy, scipy,
h5py, geopack). The numerical record is `traces.h5`; `traces.json` is the browser
transport. Serve this folder over HTTP. No backend or API keys are required.

T96 external field is added to a centered 30500 nT dipole. GSM in Earth radii.
Integration uses adaptive RK45, rtol=2e-6, atol=1e-8, max step 0.16 RE.
Background field traces terminate at Earth, 25 RE, or 85 RE arc length.
FAC legs terminate at the magnetic equator and preserve their full 3-D geometry.
Seed positions, tilt and model inputs are in the generator, not inferred from data.
The optional GSM X-Z diagnostic evaluates Jy=(dBx/dz-dBz/dx)/mu0 from the
displayed model on Y=0. It is disabled by default. The inner 2 RE are masked;
outside the displayed T96 magnetopause B is set to zero so the boundary sheet
is included in the curl. Red is +Y (dusk), blue is -Y (dawn), in nA/m^2.
The colored line circuits remain illustrative; no conductances or observed
measurements are plotted.

Included: R1/R2 FAC in both hemispheres, Pedersen connectors, partial/symmetric
ring current, and cross-tail current whose return is on the T96 magnetopause.
Chapman-Ferraro shielding and the tail return share one boundary-current layer.
R1 no longer disappears at the magnetic equator. Following Figure 4 of
Ganushkina et al. (2018), the same representative R1 sheet has two explicitly
labeled alternative magnetospheric paths: (1) a compact open-field connection
directly through the dayside magnetopause/solar-wind generator and (2) a long
closed-field route that travels antisunward in the plasma sheet, turns around
the far-tail edge, and returns earthward on the high-latitude magnetopause.
The ionospheric segment is shared by either alternative. The second path is not
a dusk-to-dawn bridge at one X location. These perpendicular paths are
schematic because tracing B does not uniquely solve them, and the visualization
does not claim that R1 simply merges into Chapman-Ferraro current. See
https://doi.org/10.1002/2017RG000590.
Included additionally: an equivalent substorm wedge in both hemispheres, with
T96-traced upward premidnight and downward postmidnight FAC, a westward auroral
electrojet at 66 degrees magnetic latitude, and schematic eastward nightside
equatorial closure. It is switchable independently of R1/R2; neither preset
predicts whether a substorm occurs. See https://doi.org/10.1007/s11214-014-0124-9.
Omitted: general Hall electrojets, NBZ/cusp and time-dependent coupling.
Chapman-Ferraro directions are evaluated on T96's pressure-scaled sigma=1.08
magnetopause using K parallel to (B_inside - B_outside) cross outward normal.
Here B_inside is the centered dipole plus T96 dipole shielding; B_outside=0 is
an idealized unmagnetized-sheath assumption. These are shielding-current
streamlines, not full model current densities. Streamlines are truncated at
X=-25 RE or the integration-length limit, not claimed to be closed circuits.
The translucent shell is the same T96 boundary, not a separate Shue model.

BCBFs are shown separately as dashed white plasma-flow channels, not
conventional-current arrows (Archer and Knudsen,
https://doi.org/10.1002/2017JA024577). Illustrative 66.3-66.7 degree channels
span 18-23 MLT westward and 1-6 MLT eastward in both hemispheres. This fixed
geometry does not predict observed widths, speeds, electric fields or onset.
The existing schematic R1/R2 footpoints remain at 70/63 degrees. BCBFs are
not identified with the substorm wedge or its electrojet.
This is a global-system schematic, not a complete quantitative magnetosphere model.

Conventional R1 current is upward dusk, downward dawn; R2 is opposite. The
partial ring path goes westward through midnight. Cross-tail current goes dawn
to dusk (+Y). Arrow speeds are illustrative.

Deploy only this folder to `/var/www/html/currents/` on SSH alias `juha-no`;
do not replace the homepage or BELLA. Validate presets, orbiting, layer toggles,
pause, and ionosphere camera on the live URL after deployment.
