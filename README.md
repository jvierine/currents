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
R2 FAC legs terminate at the magnetic equator. R1 follows T96 only through
the inner rising leg, to about 4 RE, then joins the schematic outer sheet.
Seed positions, tilt and model inputs are in the generator, not inferred from data.
The optional GSM X-Z diagnostic evaluates Jy=(dBx/dz-dBz/dx)/mu0 from the
displayed model on Y=0. It is disabled by default and spans X=-27..13 RE and
Z=-30..30 RE. The inner 5 RE are excluded because the internal-field curl is
not a magnetospheric current system. Outside the displayed T96 magnetopause B
is set to zero so the boundary sheet is included in the curl. Red is +Y
(dusk), blue is -Y (dawn), in nA/m^2.
The colored line circuits remain illustrative; no conductances or observed
measurements are plotted.

Included: R1/R2 FAC in both hemispheres, Pedersen connectors, partial/symmetric
ring current, and cross-tail current whose return is on the T96 magnetopause.
Chapman-Ferraro shielding and the tail return share one boundary-current layer.
R1 no longer disappears at the magnetic equator or self-closes across the
polar cap. Following the red Region 1 system in Figure 7 of Ganushkina et al.
(2018), the dusk and dawn R1 FACs join through a schematic high-latitude
Chapman-Ferraro/tail-boundary transition. At dawn and dusk, ionospheric
Pedersen segments connect R1 to R2. R2 joins the partial ring current, so the
R1 and R2 systems form one continuous circuit. This is the explicit topology:
R1 dusk FAC -> high-latitude boundary -> R1 dawn FAC -> dawn Pedersen -> R2
dawn FAC -> partial ring -> R2 dusk FAC -> dusk Pedersen -> R1 dusk FAC.
The perpendicular R1 branch is schematic because tracing B does not uniquely
solve generator-region current closure. See
https://doi.org/10.1002/2017RG000590.
The red ribbon comprises five nested paths extending from the rising FAC
legs to the outer high-latitude boundary behind the cusp. It never descends
to the equator and bends back up. The quiet ribbon spans outer X=-4..-13 RE;
the active ribbon expands to X=+2..-13 RE. These extents are illustrative.
The Figure 7 button isolates red R1 and green Chapman-Ferraro sheets.
The substorm wedge follows Ganushkina et al. (2018), Figure 9b, as a diversion
of the tail current through the ionosphere. Five neighboring T96 traces at
65.2–66.4 degrees form pink FAC ribbons: downward dawn/postmidnight, westward
auroral electrojet, upward dusk/premidnight. Blue tail-current feeders run
dawn to dusk and join a magnetopause return. This replaces the former isolated
reverse equatorial arc (a perturbation-loop depiction, not the figure's full
diverted circuit). Tail connections and ribbon widths are schematic, not
T96-derived current magnitudes. The wedge toggle includes its electrojet and
blue feeds/return; Figure 9 in the view controls, or `?view=figure9`, highlights
the system. `build_traces.py` generates all geometry in `traces.h5` and
`traces.json`. Neither preset predicts a substorm.
Reference: https://doi.org/10.1002/2017RG000590, section 4.3 and Figure 9.
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
The schematic R1 footpoints are at 74 degrees (quiet) and 72 (active),
with R2 at 63 degrees. BCBFs are
not identified with the substorm wedge or its electrojet.
This is a global-system schematic, not a complete quantitative magnetosphere model.

Conventional R1 current is upward dusk, downward dawn; R2 is opposite. The
partial ring path goes westward through midnight. Cross-tail current goes dawn
to dusk (+Y). Arrow speeds are illustrative.

Deploy only this folder to `/var/www/html/currents/` on SSH alias `juha-no`;
do not replace the homepage or BELLA. Validate presets, orbiting, layer toggles,
pause, and ionosphere camera on the live URL after deployment.
