# ======================================================================
# HARWICH-FORMATION BOREHOLES – INTERACTIVE VIEWER
# ======================================================================
####

import streamlit as st
import pandas as pd
import numpy as np
import os
import pydeck as pdk
from pyproj import Transformer
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import random

# ----------------------------------------------------------------------
# 0. Streamlit page
# ----------------------------------------------------------------------
st.set_page_config(page_title="Harwich Formation Boreholes", layout="wide")
st.title("Boreholes that Intercept the Harwich Formation")

# ----------------------------------------------------------------------
# 1. Pick an AGS workbook in the current directory
# ----------------------------------------------------------------------
xlsx_files = [f for f in os.listdir()
              if f.lower().endswith(".xlsx") and not f.startswith("~$")]

excel_file = st.selectbox(
    "Available .xlsx files in this folder:",
    options=xlsx_files,
    index=0 if xlsx_files else None,
    placeholder="Choose an AGS workbook"
)

if not excel_file:
    st.warning("No AGS Excel (.xlsx) files found in this folder.")
    st.stop()

# ----------------------------------------------------------------------
# 2. Read LOCA & GEOL sheets
# ----------------------------------------------------------------------
LOCA_SHEET = "LOCA - AGS"
GEOL_SHEET = "GEOL - AGS"
GEOL_COLS  = ["LOCA_ID",
              "LOCA_GL",
              "GEOL_TOP",
              "GEOL_BASE",
              "COWI GEOL_2",
              "COWI GEOL_4"]

try:
    # ---- LOCA ----
    loca_df = pd.read_excel(excel_file, sheet_name=LOCA_SHEET, header=3)
    loca_df.columns = loca_df.columns.str.strip()
    loca_df = loca_df[["LOCA_ID", "LOCA_NATE", "LOCA_NATN"]]

    loca_df[["LOCA_NATE", "LOCA_NATN"]] = loca_df[["LOCA_NATE", "LOCA_NATN"]].apply(
        pd.to_numeric, errors="coerce")
    loca_df = loca_df.dropna(subset=["LOCA_NATE", "LOCA_NATN"])

    transformer = Transformer.from_crs("epsg:27700", "epsg:4326", always_xy=True)

    def osgb_to_latlon(row):
        try:
            lon, lat = transformer.transform(float(row["LOCA_NATE"]),
                                             float(row["LOCA_NATN"]))
            return pd.Series({"LAT": lat, "LON": lon})
        except Exception:
            return pd.Series({"LAT": np.nan, "LON": np.nan})

    loca_df[["LAT", "LON"]] = loca_df.apply(osgb_to_latlon, axis=1)
    loca_df = loca_df.dropna(subset=["LAT", "LON"])

    # ---- GEOL ----
    geol_df = pd.read_excel(excel_file,
                            sheet_name=GEOL_SHEET,
                            usecols=GEOL_COLS)
    geol_df.columns = geol_df.columns.str.strip()

except Exception as e:
    st.error(f"Error loading AGS file: {e}\n\n"
             "Check sheet names, header rows and column names.")
    st.stop()

# ----------------------------------------------------------------------
# 3. Keep ONLY boreholes that mention “Harwich” in GEOL_2 or GEOL_4
# ----------------------------------------------------------------------
def contains_harwich(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().str.contains("HARWICH", na=False)

harwich_mask = (
    contains_harwich(geol_df["COWI GEOL_2"]) |
    contains_harwich(geol_df["COWI GEOL_4"])
)
harwich_ids = geol_df.loc[harwich_mask, "LOCA_ID"].unique()

if len(harwich_ids) == 0:
    st.warning(
        "No intervals contain the word 'Harwich'.\n\n"
        "Distinct values found in COWI GEOL_2 and GEOL_4 are listed below."
    )
    col1, col2 = st.columns(2)
    with col1:
        st.write("COWI GEOL_2")
        st.dataframe(pd.DataFrame(geol_df["COWI GEOL_2"].unique(),
                                  columns=["COWI GEOL_2"]),
                     hide_index=True)
    with col2:
        st.write("COWI GEOL_4")
        st.dataframe(pd.DataFrame(geol_df["COWI GEOL_4"].unique(),
                                  columns=["COWI GEOL_4"]),
                     hide_index=True)
    st.stop()

loca_df = loca_df[loca_df["LOCA_ID"].isin(harwich_ids)].reset_index(drop=True)
geol_df = geol_df[geol_df["LOCA_ID"].isin(harwich_ids)]

# ----------------------------------------------------------------------
# 4. Colour palette for lithology (COWI GEOL_4)
# ----------------------------------------------------------------------
unique_desc = geol_df["COWI GEOL_4"].dropna().unique()
palette = (list(mcolors.TABLEAU_COLORS.values()) +
           list(mcolors.XKCD_COLORS.values()))
random.shuffle(palette)
DESC_COLOUR = {d: palette[i % len(palette)] for i, d in enumerate(unique_desc)}
DEFAULT_COLOUR = "#cccccc"

# ----------------------------------------------------------------------
# 5. Layout – right: map, left: lithology stick
# ----------------------------------------------------------------------
left_col, right_col = st.columns([1, 2])

# --------------------------- Map --------------------------------------
with right_col:
    st.subheader("Harwich borehole locations")

    view_state = pdk.ViewState(
        longitude=loca_df["LON"].mean(),
        latitude=loca_df["LAT"].mean(),
        zoom=13,
        pitch=0,
    )

    base_layer = pdk.Layer(
        "ScatterplotLayer",
        data=loca_df,
        pickable=False,
        get_position="[LON, LAT]",
        get_radius=8,
        get_fill_color="[200, 30, 0, 140]",
        get_line_color="[0, 0, 0]",
        line_width_min_pixels=1,
    )

# ---------------------- Lithology stick & controls --------------------
with left_col:
    st.subheader("Lithology")

    # ---- dropdown + manual input -------------------------------------
    id_options = sorted(loca_df["LOCA_ID"].astype(str).tolist())
    dropdown_id = st.selectbox(
        "Choose a borehole (Harwich only):",
        [""] + id_options,
        index=0,
    )
    manual_id = st.text_input("…or type a borehole ID:")
    chosen_id = (manual_id.strip() or dropdown_id).strip()

    # default – nothing selected
    selected_row = pd.DataFrame()
    layers = [base_layer]

    if chosen_id == "":
        st.markdown("_Select or type a borehole ID to see details._")

    elif chosen_id not in id_options:
        st.error(f"Borehole '{chosen_id}' does not exist in the filtered list.")

    else:
        # ------------- build lithology plot ---------------------------
        intervals = geol_df[geol_df["LOCA_ID"] == chosen_id]\
            .sort_values("GEOL_TOP")

        numeric_cols = ["LOCA_GL", "GEOL_TOP", "GEOL_BASE"]
        intervals[numeric_cols] = intervals[numeric_cols].apply(
            pd.to_numeric, errors="coerce")
        intervals = intervals[np.isfinite(intervals[numeric_cols]).all(axis=1)]

        if intervals.empty:
            st.warning("No plottable intervals for this borehole.")
        else:
            # axis limits only around logged part
            ground_level = intervals.iloc[0]["LOCA_GL"]
            top_od = ground_level - intervals["GEOL_TOP"]
            base_od = ground_level - intervals["GEOL_BASE"]

            axis_top = top_od.max() + 0.5
            axis_bottom = base_od.min() - 0.5

            fig, ax = plt.subplots(figsize=(2.2, 6))

            for y_top, y_base, desc in zip(top_od, base_od,
                                           intervals["COWI GEOL_4"]):
                colour = DESC_COLOUR.get(desc, DEFAULT_COLOUR)
                thickness = y_top - y_base
                ax.add_patch(
                    patches.Rectangle(
                        (0.3, y_base), 0.4, thickness,
                        facecolor=colour, edgecolor="black", linewidth=0.6
                    )
                )

                txt = str(desc)
                txt = txt[:12] + "…" if len(txt) > 12 else txt
                if thickness >= 1.0:
                    ax.text(0.5, (y_top + y_base) / 2, txt,
                            ha="center", va="center", fontsize=8)
                else:
                    ax.text(0.75, (y_top + y_base) / 2, txt,
                            ha="left", va="center", fontsize=7)

            # dynamic tick spacing
            span = axis_top - axis_bottom
            tick_step = 1 if span <= 5 else (2 if span <= 20 else 5)
            ticks = np.arange(np.floor(axis_bottom / tick_step) * tick_step,
                              axis_top + tick_step, tick_step)

            ax.set_ylim(axis_top, axis_bottom)   # ground level on top
            ax.set_xlim(0, 1.1)
            ax.set_yticks(ticks)
            ax.set_ylabel("Elevation (m OD)")
            ax.invert_yaxis()

            ax.xaxis.set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)

            st.pyplot(fig)

            # interval table
            st.markdown("**Interval table**")
            st.dataframe(
                pd.DataFrame({
                    "Top Elev (m OD)": top_od,
                    "Base Elev (m OD)": base_od,
                    "Description": intervals["COWI GEOL_4"]
                }),
                hide_index=True,
            )

            # highlight on map
            selected_row = loca_df[loca_df["LOCA_ID"] == chosen_id]
            highlight_layer = pdk.Layer(
                "ScatterplotLayer",
                data=selected_row,
                get_position="[LON, LAT]",
                get_radius=15,
                get_fill_color="[0, 0, 0, 255]",
                get_line_color="[255, 255, 255]",
                line_width_min_pixels=2,
                pickable=False,
            )
            layers = [base_layer, highlight_layer]

# ----------------------- Render the map --------------------------------
with right_col:
    deck = pdk.Deck(
        map_style="mapbox://styles/mapbox/light-v9",
        initial_view_state=view_state,
        layers=layers,
        tooltip={"text": "LOCA_ID: {LOCA_ID}"},
    )
    st.pydeck_chart(deck, use_container_width=True)