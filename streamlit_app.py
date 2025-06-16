import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

# --- PAGE CONFIG --- #
st.set_page_config(
    page_title="Hardening Soil Model",
    layout="wide"
)

# --- SIDEBAR --- #
sb = st.sidebar
max_depth = sb.slider(
    "Maximum depth for figure (m)", min_value=10, max_value=100, value=20
)
sb.markdown(
    "Read about Hardening Soil in [this blog post.](https://berkdemir.github.io/2021/05/11/Hardening-Soil-Model/)"
)
sb.caption(
    "Boussinesq only for vertical stress; horizontal = K0 * vertical."
)

# --- MAIN TITLE --- #
st.title("Hardening Soil Model | Effect of Power m")

# --- PARAMETER INPUTS --- #
columns = st.columns([1, 1, 3])
with columns[0]:
    B = st.slider("Width of Foundation (m)", 1, 100, 20)
    L = st.slider("Length of Foundation (m)", 1, 100, 40)
    q = st.slider("Pressure on Foundation (kPa)", 0, 500, 100)
    UW = st.slider("Unit Weight of Soil (kN/m3)", 0, 25, 20)
    K0 = st.slider("Coefficient of Lateral Earth Pressure", 0.0, 2.0, 0.5)

with columns[1]:
    m = st.slider("Power m", 0.0, 1.0, 0.55)
    pref = st.slider("Reference Pressure, pref (kPa)", 10, 500, 100)
    E50ref = st.slider("E50ref (MPa)", 10, 500, 50)
    c = st.slider("Cohesion (kPa)", 0, 200, 20)
    phi = st.slider("Angle of Friction (deg)", 0, 45, 30)

# --- CORE CALCULATIONS --- #
def boussinesq(L, B, z, q):
    mL = L / B
    b = B / 2
    n = z / b
    I = (
        (
            mL * n * (1 + mL**2 + 2*n**2) /
            np.sqrt(1 + mL**2 + n**2) /
            (1 + n**2) /
            (mL**2 + n**2)
            + np.arcsin(mL / (np.sqrt(mL**2 + n**2) * np.sqrt(1 + n**2)))
        ) * 2 / np.pi
    )
    return q * I

def Ecalc(E50ref, q, K0, m, pref, c, phi):
    return E50ref * (
        (c * np.cos(np.radians(phi)) + q * K0 * np.sin(np.radians(phi))) /
        (c * np.cos(np.radians(phi)) + pref * np.sin(np.radians(phi)))
    ) ** m

depth = np.arange(0, max_depth+1, 1)
qred = []
E50_load = []
E50_load_m0 = []
E50_load_m05 = []
E50_load_m1 = []

for i in depth:
    q_eff = boussinesq(L, B, i, q) + UW * i
    qred.append(q_eff)
    E50_load.append(Ecalc(E50ref, q_eff, K0, m, pref, c, phi))
    E50_load_m0.append(Ecalc(E50ref, q_eff, K0, 0, pref, c, phi))
    E50_load_m05.append(Ecalc(E50ref, q_eff, K0, 0.5, pref, c, phi))
    E50_load_m1.append(Ecalc(E50ref, q_eff, K0, 1, pref, c, phi))

# --- PLOTTING --- #
fig, ax = plt.subplots(1, 2, figsize=(10, 6), sharey=True)
ax[0].plot(qred, depth, label="Pressure with Depth", color="red")
ax[0].set_xlabel("Pressure (kPa)")
ax[0].set_ylabel("Depth (m)")
ax[0].invert_yaxis()
ax[0].legend(fontsize=8)
ax[0].set_title("Pressure Difference with Depth (Boussinesq)")

ax[1].plot(E50_load, depth, label=f"m = {m:.2f}", color="red")
ax[1].plot(E50_load_m0, depth, label="m = 0", color="red", ls="--", lw=1)
ax[1].plot(E50_load_m05, depth, label="m = 0.5", color="green", ls="--", lw=0.5)
ax[1].plot(E50_load_m1, depth, label="m = 1", color="gray", ls="--", lw=1)
ax[1].set_xlabel("E50 (MPa)")
ax[1].legend(fontsize=8)
ax[1].set_title("E50 with Depth")
ax[1].fill_betweenx(depth, E50_load_m1, E50_load_m0, color="gray", alpha=0.1)

fig.suptitle("Effect of Load and Power m on E50", fontsize=16)
fig.tight_layout()

with columns[2]:
    st.pyplot(fig)