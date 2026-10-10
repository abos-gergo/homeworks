const { vec3 } = require('gl-matrix');

// Adatok
const l1 = 0.06;
const l2 = 0.15;
const r3 = 0.04;
const R = 3 * r3; // 0.12
const omega1z = 5.0;
const phi_deg = 50.0;
const phi3_deg = 60.0;

const deg2rad = (deg) => deg * Math.PI / 180.0;
const rad2deg = (rad) => rad * 180.0 / Math.PI;

const phi = deg2rad(phi_deg);
const phi3 = deg2rad(phi3_deg);

// Pozíciók (vektorok)
const rA = vec3.fromValues(0, 0, 0);
const rB = vec3.fromValues(l1 * Math.cos(phi), l1 * Math.sin(phi), 0);

// A korábbi számításod alapján a (Dx, Dy) a talppontot (érintkezési pontot) adta meg
const Dy = rB[1];
const Cy = Dy + r3 + r3 * Math.cos(phi3);
const Cx = rB[0] + Math.sqrt(l2 ** 2 - (Cy - rB[1]) ** 2);
const Dx = Cx + r3 * Math.sin(phi3);

const rC = vec3.fromValues(Cx, Cy, 0);
const rD = vec3.fromValues(Dx, Dy, 0); // K (érintkezési pont a talajon)
const rS3 = vec3.fromValues(Dx, Dy + r3, 0); // 3-as test súlypontja (kerék közepe)
const rS2 = vec3.fromValues((rB[0] + rC[0]) / 2, (rB[1] + rC[1]) / 2, 0); // 2-es rúd súlypontja

// --- 2. FELADAT: omega2 és vS2 ---
// Pólus meghatározása (P2)
const bx = rB[0], by = rB[1];
const dx = rD[0], dy = rD[1];
const cx = rC[0], cy = rC[1];
const det = by * (dx - cx) - bx * (dy - cy);
const crossDC = dx * cy - dy * cx;
const Px = (bx * crossDC) / det;
const Py = (by * crossDC) / det;
const rP2 = vec3.fromValues(Px, Py, 0);

const omega1 = vec3.fromValues(0, 0, omega1z);
const vB = vec3.create();
vec3.cross(vB, omega1, rB);

const rP2B = vec3.create(); vec3.subtract(rP2B, rB, rP2);
const cross_rP2B_vB = vec3.create(); vec3.cross(cross_rP2B_vB, rP2B, vB);
const omega2z = cross_rP2B_vB[2] / (rP2B[0] ** 2 + rP2B[1] ** 2);
const omega2 = vec3.fromValues(0, 0, omega2z);

const rP2S2 = vec3.create(); vec3.subtract(rP2S2, rS2, rP2);
const vS2 = vec3.create(); vec3.cross(vS2, omega2, rP2S2);
const vS2_mag = vec3.length(vS2);

// vC és omega3 számítása (a webes űrlaphoz)
const rP2C = vec3.create(); vec3.subtract(rP2C, rC, rP2);
const vC = vec3.create(); vec3.cross(vC, omega2, rP2C);
const omega3z = vC[1] / (rC[0] - rS3[0]);

// --- 4. FELADAT: eps2 és aS2 (valamint eps3) ---
const aB = vec3.fromValues(-(omega1z ** 2) * rB[0], -(omega1z ** 2) * rB[1], 0);

const vS3x = -omega3z * r3;
const aS3y = -(vS3x ** 2) / (R + r3); // A gördülő kerék középpontjának normális gyorsulása

const rBCx = rC[0] - rB[0], rBCy = rC[1] - rB[1];
const rS3Cx = rC[0] - rS3[0], rS3Cy = rC[1] - rS3[1];

// A kinematikai kényszerekből adódó lineáris egyenletrendszer
const A11 = -rBCy; const A12 = r3 + rS3Cy;
const B1_ = -aB[0] + (omega2z ** 2) * rBCx - (omega3z ** 2) * rS3Cx;
const A21 = rBCx; const A22 = -rS3Cx;
const B2_ = aS3y - aB[1] + (omega2z ** 2) * rBCy - (omega3z ** 2) * rS3Cy;

const detA = A11 * A22 - A12 * A21;
const eps2z = (B1_ * A22 - B2_ * A12) / detA;
const eps3z = (A11 * B2_ - A21 * B1_) / detA;
const eps2 = vec3.fromValues(0, 0, eps2z);

const rBS2 = vec3.create(); vec3.subtract(rBS2, rS2, rB);
const aS2 = vec3.create(); vec3.cross(aS2, eps2, rBS2);
aS2[0] += aB[0] - (omega2z ** 2) * rBS2[0];
aS2[1] += aB[1] - (omega2z ** 2) * rBS2[1];
const aS2_mag = vec3.length(aS2);

// aC számítás rajzoláshoz
const aC = vec3.create(); vec3.cross(aC, eps2, vec3.fromValues(rBCx, rBCy, 0));
aC[0] += aB[0] - (omega2z ** 2) * rBCx;
aC[1] += aB[1] - (omega2z ** 2) * rBCy;

// --- 6. FELADAT: gyorsulásszög és gyorsuláspólus ---
const alpha2_rad = Math.atan2(Math.abs(eps2z), omega2z ** 2);
const alpha2_deg = rad2deg(alpha2_rad);

const eps2_cross_aB = vec3.create(); vec3.cross(eps2_cross_aB, eps2, aB);
const denom = (omega2z ** 4) + (eps2z ** 2);
const rBG2 = vec3.fromValues(
    ((omega2z ** 2) * aB[0] + eps2_cross_aB[0]) / denom,
    ((omega2z ** 2) * aB[1] + eps2_cross_aB[1]) / denom,
    0
);
const rG2 = vec3.create(); vec3.add(rG2, rB, rBG2);

// --- 7. FELADAT: aS2t és aS2n ---
const e_t = vec3.fromValues(vS2[0] / vS2_mag, vS2[1] / vS2_mag, 0);
const aS2t_scalar = vec3.dot(aS2, e_t);
const aS2t_mag = Math.abs(aS2t_scalar);
const aS2t = vec3.fromValues(aS2t_scalar * e_t[0], aS2t_scalar * e_t[1], 0);

const aS2n = vec3.create(); vec3.subtract(aS2n, aS2, aS2t);
const aS2n_mag = vec3.length(aS2n);

// --- 8. FELADAT: rho_S2 ---
const rho_S2 = (vS2_mag ** 2) / aS2n_mag;

// Konzol output generálása
console.log("=========================================");
console.log("          HW CHK WEB FORM DATA           ");
console.log("=========================================");
console.log(`omega_2z : ${omega2z.toFixed(5)}  [rad/s]`);
console.log(`eps_2z   : ${eps2z.toFixed(5)}  [rad/s^2]`);
console.log(`omega_3z : ${omega3z.toFixed(5)}  [rad/s]`);
console.log(`eps_3z   : ${eps3z.toFixed(5)}  [rad/s^2]`);
console.log(`v_S2     : ${vS2_mag.toFixed(5)}  [m/s]`);
console.log(`a_S2     : ${aS2_mag.toFixed(5)}  [m/s^2]`);
console.log(`|alpha_2|: ${alpha2_deg.toFixed(5)}  [fok]`);
console.log(`|a_S2t|  : ${aS2t_mag.toFixed(5)}  [m/s^2]`);
console.log(`|a_S2n|  : ${aS2n_mag.toFixed(5)}  [m/s^2]`);
console.log(`rho_S2   : ${rho_S2.toFixed(5)}  [m]`);
console.log("=========================================\n");

console.log("--- Vektorok a rajzolós feladatokhoz (1, 3, 5. feladat) ---");
console.log(`[Pozíciók]`);
console.log(`A  = (${rA[0].toFixed(5)}, ${rA[1].toFixed(5)}) m`);
console.log(`B  = (${rB[0].toFixed(5)}, ${rB[1].toFixed(5)}) m`);
console.log(`C  = (${rC[0].toFixed(5)}, ${rC[1].toFixed(5)}) m`);
console.log(`K (érintők) = (${rD[0].toFixed(5)}, ${rD[1].toFixed(5)}) m`);
console.log(`S3 = (${rS3[0].toFixed(5)}, ${rS3[1].toFixed(5)}) m`);
console.log(`S2 = (${rS2[0].toFixed(5)}, ${rS2[1].toFixed(5)}) m`);

console.log(`\n[Sebességpólus]`);
console.log(`P2 = (${rP2[0].toFixed(5)}, ${rP2[1].toFixed(5)}) m`);

console.log(`\n[Sebességek]`);
console.log(`vB  = (${vB[0].toFixed(5)}, ${vB[1].toFixed(5)}) m/s`);
console.log(`vC  = (${vC[0].toFixed(5)}, ${vC[1].toFixed(5)}) m/s`);
console.log(`vS2 = (${vS2[0].toFixed(5)}, ${vS2[1].toFixed(5)}) m/s`);

console.log(`\n[Gyorsuláspólus]`);
console.log(`G2 = (${rG2[0].toFixed(5)}, ${rG2[1].toFixed(5)}) m`);

console.log(`\n[Gyorsulások]`);
console.log(`aB  = (${aB[0].toFixed(5)}, ${aB[1].toFixed(5)}) m/s^2`);
console.log(`aC  = (${aC[0].toFixed(5)}, ${aC[1].toFixed(5)}) m/s^2`);
console.log(`aS2 = (${aS2[0].toFixed(5)}, ${aS2[1].toFixed(5)}) m/s^2`);
console.log(`aS2t (vektor) = (${aS2t[0].toFixed(5)}, ${aS2t[1].toFixed(5)}) m/s^2`);
console.log(`aS2n (vektor) = (${aS2n[0].toFixed(5)}, ${aS2n[1].toFixed(5)}) m/s^2`);
