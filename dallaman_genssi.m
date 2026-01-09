function model = dallaman_genssi()
% GenSSI structural identifiability model
% Exact Dalla Man 2007 equations + CGM sensor
% MATLAB-safe symbolic names

%% =========================
% STATES
% =========================
syms Gp Gt Il Ip Qsto1 Qsto2 Qgut I1 Id X Ipo Y real
model.sym.x = [Gp; Gt; Il; Ip; Qsto1; Qsto2; Qgut; I1; Id; X; Ipo; Y];
model.sym.Name = 'dallaman_genssi';   % REQUIRED FOR GenSSI REPORTING
model.sym.Nder = 4;                   % number of Lie derivatives
% =========================
% INITIAL CONDITIONS
% =========================
syms Gp0 Gt0 Il0 Ip0 Qsto10 Qsto20 Qgut0 I10 Id0 X0 Ipo0 Y0 real
model.sym.x0 = [ ...
    Gp0; Gt0; Il0; Ip0; ...
    Qsto10; Qsto20; Qgut0; ...
    I10; Id0; X0; Ipo0; Y0 ];

%% =========================
% PARAMETERS
% =========================
syms V_G k_1 k_2 G_b V_I ...
     m_1 m_2 m_4 m_5 m_6 ...
     k_max k_min k_abs k_gri ...
     f b d BW ...
     k_p1 k_p2 k_p3 k_p4 ...
     k_i U_ii ...
     V_m0 V_mX K_m0 ...
     p_2U part_frac K ...
     alpha_y beta_y gamma_i ...
     k_e1 k_e2 D ...
     tau_cgm real

model.sym.p = [ ...
    V_G k_1 k_2 G_b V_I ...
    m_1 m_2 m_4 m_5 m_6 ...
    k_max k_min k_abs k_gri ...
    f b d BW ...
    k_p1 k_p2 k_p3 k_p4 ...
    k_i U_ii ...
    V_m0 V_mX K_m0 ...
    p_2U part_frac K ...
    alpha_y beta_y gamma_i ...
    k_e1 k_e2 D ...
    tau_cgm ]';
model.sym.Par = model.sym.p;   % REQUIRED for GenSSI report

%% =========================
% INPUT
% =========================
syms meal real
model.sym.u = meal;

%% =========================
% AUXILIARY EXPRESSIONS
% =========================
EGP = k_p1 - k_p2*Gp - k_p3*Id - k_p4*Ipo;

V_mmax = (1 - part_frac) * (V_m0 + V_mX*X);
Uidm = V_mmax * Gt / (K_m0 + Gt);

G = Gp / V_G;
I = Ip / V_I;

HE = -m_5 * (gamma_i*Ipo) + m_6;
m_3 = m_1*HE/(1-HE);

Qsto = Qsto1 + Qsto2;
Ra = f * k_abs * Qgut / BW;

aa = 5/(2*(1-b)*D);
cc = 5/(2*d*D);

k_empt = k_min + (k_max - k_min)/2 * ...
    (tanh(aa*(Qsto - b*D)) - tanh(cc*(Qsto - d*D)) + 2);

Spo = Y + K*(EGP + Ra - U_ii - k_1*Gp + k_2*Gt)/V_G;

%% =========================
% DIFFERENTIAL EQUATIONS
% =========================
dGp    = EGP + Ra - U_ii - k_1*Gp + k_2*Gt;
dGt    = -Uidm + k_1*Gp - k_2*Gt;
dIl    = -m_1*Il - m_3*Il + m_2*Ip + gamma_i*Ipo;
dIp    = -m_2*Ip - m_4*Ip + m_1*Il;

dQsto1 = -k_gri*Qsto1;
dQsto2 = -k_empt*Qsto2 + k_gri*Qsto1;
dQgut  = -k_abs*Qgut + k_empt*Qsto2;

dI1    = -k_i*(I1 - I);
dId    = -k_i*(Id - I1);

dX     = -p_2U*X + p_2U*(I - (Ip/V_I));
dIpo   = -gamma_i*Ipo + Spo;
dY     = -alpha_y*(Y - beta_y*(G - G_b));

model.sym.f = [ ...
    dGp; dGt; dIl; dIp;
    dQsto1; dQsto2; dQgut;
    dI1; dId; dX; dIpo; dY ];

model.sym.xdot = model.sym.f;   % GenSSI-required alias


%% =========================
% OUTPUT (CGM SENSOR)
% =========================
model.sym.h = Y;
model.sym.y = model.sym.h;   % GenSSI-required alias


end
