import numpy as np
from shared.LMCS import glob_util
from shared.paths import load_paths
import os
import pandas as pd


REGIONS = glob_util.REGIONS
style = glob_util.REGION_STYLE


paths = load_paths()
ROOT = paths.datasets["obs_mcs_5000km2_tables_ERA"]
print("Using ROOT:", ROOT)
YEARS = range(2000, 2020)

def make_df():

    def read_meta(region):
        files = [
            f"{ROOT}/storm_metadata/region={region}/year={y}/part.parquet"
            for y in YEARS
        ]

        files = [f for f in files if os.path.isfile(f)]

        return pd.concat(
            [pd.read_parquet(f) for f in files],
            ignore_index=True
        )


    def read_var(var, region):
        files = [
            f"{ROOT}/{var}/region={region}/year={y}/part.parquet"
            for y in YEARS
        ]

        files = [f for f in files if os.path.isfile(f)]

        return pd.concat(
            [
                pd.read_parquet(
                    f,
                    columns=["storm_id", "rel_hour", "value"]
                )
                for f in files
            ],
            ignore_index=True
        )


    def build_all_regions(variables):

        dfs = []

        for region in REGIONS:

            meta = read_meta(region)

            if meta is None:
                continue

            df = None

            for var in variables:

                tmp = read_var(var, region)

                tmp = tmp.rename(columns={"value": var})

                if df is None:
                    df = tmp
                else:
                    df = df.merge(
                        tmp[["storm_id", "rel_hour", var]],
                        on=["storm_id", "rel_hour"],
                        how="outer",
                    )

            df = df.merge(meta, on="storm_id", how="left")

            df["region"] = region

            dfs.append(df)

        return pd.concat(dfs, ignore_index=True)


    def add_rh_from_q_t(
        df,
        q_col="q650",
        t_col="t650",
        rh_col="rh650",
        p_hpa=650.0,
    ):
        """
        Add relative humidity (%) from specific humidity q and temperature T.

        q must be specific humidity in kg/kg.
        T must be temperature in K.
        p_hpa is pressure in hPa.
        """

        q = df[q_col].astype(float)
        T = df[t_col].astype(float)

        # Convert pressure to Pa
        p = p_hpa * 100.0

        # Specific humidity q -> vapour pressure e
        # q = epsilon * e / (p - (1 - epsilon) * e)
        # rearranged:
        epsilon = 0.622
        e = q * p / (epsilon + (1.0 - epsilon) * q)

        # Saturation vapour pressure over liquid water, Bolton-style formula
        # T in K, output in Pa
        T_c = T - 273.15
        es = 611.2 * np.exp((17.67 * T_c) / (T_c + 243.5))

        rh = 100.0 * e / es

        df[rh_col] = rh

        return df




    def add_low_level_moisture_metrics(df):
        """
        Adds low/mid-level moisture diagnostics from q925, q850, q750.

        Assumes q is specific humidity in kg/kg.
        Output IWV is kg/m2, equivalent to mm.
        """

        g = 9.80665

        q925 = df["q925"].astype("float32")
        q850 = df["q850"].astype("float32")
        q750 = df["q750"].astype("float32")

        # pressure-layer thicknesses in Pa
        dp_750_850 = (850 - 750) * 100.0
        dp_850_925 = (925 - 850) * 100.0
        dp_total = (925 - 750) * 100.0

        iwv = (
            0.5 * (q750 + q850) * dp_750_850
            + 0.5 * (q850 + q925) * dp_850_925
        ) / g

        df["iwv_925_750"] = iwv.astype("float32")

        # q-like pressure-weighted mean over the layer
        df["qmean_925_750"] = (df["iwv_925_750"] * g / dp_total).astype("float32")

        return df


        ###########################################

    df_all = build_all_regions(["tcwv", "shear100_650", "u650", "ushear100_650","q925",
                                "q850","q750","q650","t650", "cape"])
    
    df_all = df_all[((df_all["ushear100_650"] > 0) & (df_all["u650"]>0)) | ((df_all["ushear100_650"] < 0)  & (df_all["u650"]< 0))
                    & (df_all["tcwv"] > 10) & (df_all["precipitation_max"] >= 3)]


    df_all = add_rh_from_q_t(df_all,q_col="q650",t_col="t650",
        rh_col="rh650",p_hpa=650.0)

    df_all["rh650"] = df_all["rh650"].clip(0, 100)


    df_all = add_low_level_moisture_metrics(df_all)
    df_all = df_all.drop(columns=["q850", "q750", "t650", "q650"], errors="ignore")

    #df_all = build_all_regions(["tcwv", "shear100_650", "u650", "ushear100_650","q925"])

    out_file = paths.shear_data + "/reduced_shear_study_data.parquet"

    df_all.to_parquet(out_file, index=False)
    print(f"Saved combined DataFrame to {out_file}")