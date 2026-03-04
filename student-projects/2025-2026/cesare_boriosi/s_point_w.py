import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
import matplotlib.pyplot as plt
# from scipy.optimize import differential_evolution

# ---------------------------
# Utility metric functions
# ---------------------------
def mape_eps(y_true, y_pred, eps=1e-8):
    """MAPE with |y| in denominator and additive eps:
       MAPE = 100/n * sum |(y - yhat)| / (|y| + eps)
    """
    denom = np.abs(y_true) + eps
    return 100.0 * np.mean(np.abs((y_true - y_pred) / denom))

def mae(y_true, y_pred):
    return np.mean(np.abs(y_true - y_pred))

def rmsd(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))

def std_residuals(y_true, y_pred):
    resid = y_true - y_pred
    return np.std(resid, ddof=1)

# importa dataframe da pickle
df_in = pd.read_pickle("Full_PBE.pkl")
print("\n-----------------------------------------------------")
print("TRAINING FOR SINGLE POINT EXTRAPOLATION")
print("-------------------------------------------------------")
print("\nFull DataFrame size:", df_in.shape)

#print(df_in.head(1))
#Full_PBE= pd.read_pickle("./Full_PBE.pkl") 
#Full_PBE

#----------------------------------------------------
#copia del df perché verrà moltiplicato
#----------------------------------------------------

df = df_in.copy()
#print("Check dataFrame size:", df.shape)



# le sole colonne che mi interessano
cols_features = [
    'Te_DZ',
    'V_eN_DZ',
    'TEE_DZ',
    ]
col_target = "Ee_QZ"
cols_to_f = [*cols_features, col_target]

df["E_res_DZ"] = df["E_DZ"] - df[cols_features].sum(axis=1)
df["Ee_QZ"] = df["E_QZ"] - df["E_res_DZ"]

# etichetta REAZIONI E ENERGIE ASSOLUTE
df["Rxn"] = df["System"].str.contains("Rxn", na=False)


# CHECK come funzina solo su una delle due metriche
df_rxn = df.loc[df["Rxn"]].copy()
df_abs = df.loc[~df["Rxn"]].copy()
print(f"Righe Energie Reazione (df_rxn): {len(df_rxn)}")
print(f"Righe Energie Assolute (df_abs): {len(df_abs)}\n\n")

print("Check dataFrame size:", df.shape)
#print(df_rxn)
# CHECK PER VEDERE STATISTICHE E PLOT DEI DATI VERGINI
check = False

cols_check = [*cols_features ,
            col_target ,
            "E_DZ" ,
            "E_QZ"]

#----------------------------------------------------
# METRICHE INIZIALI E DI ML INDIPENDENTI 
# METRICHE INIZIALI DI DIFFERENZA TRA E_DZ ED E_QZ
#----------------------------------------------------

rmsd_value_abs = rmsd(df_abs["E_QZ"], df_abs["E_DZ"])
mape_abs = mape_eps(df_abs["E_QZ"], df_abs["E_DZ"], eps = 1e-8)
rmsd_value_rxn = rmsd(df_rxn["E_QZ"], df_rxn["E_DZ"])
mape_rxn = mape_eps(df_rxn["E_QZ"], df_rxn["E_DZ"], eps = 1e-8)

# TEST DI COME FUNZIONEREBBE IL MODELLO SOLO SU ABS O SU RXN     

for label, df_test in [ ("ABSOLUTE", df_abs), ("REACTION", df_rxn)]:

    df_copy = df_test.copy().reset_index(drop=True)
    x_train = df_copy[cols_features]
    y_train = df_copy[col_target]
    model = LinearRegression(fit_intercept=False)
    model.fit(x_train, y_train)

    x_test = df_test[cols_features]
    y_pred = model.predict(x_test) + df_test["E_res_DZ"]
    y_true = df_test["E_QZ"]

    if label == ("ABSOLUTE"):
        best_rmsd_value_abs = rmsd(y_true, y_pred)
        best_mape_value_abs = mape_eps(y_true, y_pred, eps=1e-8)
    elif label == ("REACTION"):
        best_rmsd_value_rxn = rmsd(y_true, y_pred)
        best_mape_value_rxn = mape_eps(y_true, y_pred, eps=1e-8)

#----------------------------------------------------
# GRID SEARCH
#----------------------------------------------------

verbose = False

features = cols_features
target = col_target

print(f"Features usate: {features}")
print(f"Target: {target}")

print("\n  METRICS  SENZA MACHINE LAERNING DI ABSOLUTE")
print("RMSD (abs): ", rmsd_value_abs)
print("MAPE (abs): ",mape_abs)
print("\n  METRICS  SENZA MACHINE LAERNING DI RELATIVE")
print("RMSD (rxn): ", rmsd_value_rxn)
print("MAPE (rxn): ",mape_rxn)


# AGGIUNGO DEI VETTORI PER POTER FARE I GRAFICI
plot_factor = []
plot_rmsd_abs = []
plot_mape_abs = []
plot_rmsd_rxn = []
plot_mape_rxn = []

g_start = 1 
g_final = 100
g_step = 25
# modificare il numero dentro range per decidere il numero di giri
for g in range(2):
    #print("\n", g+1,  "° GIRO DI GRID SEARCH")

    if g == 0 :
        init = g_start
        fin = g_final
    else:
        steppy = (fin-init)/(g_step-1)
        init = best_factor - steppy
        fin = best_factor + steppy

    test_range = np.linspace( init, fin, g_step)
    #print("--------------------------------------------------------")
    #print(" FACTOR FROM ", test_range[0], " TO ",test_range[g_step-1])
    #test_range = [8.447368421052632] # per CARLOS
    best_deltas = np.inf
    best_factor = None
    best_betas = None

    for factor in test_range:
        
        plot_factor.append(factor)
        # inserire il dataframe da usare: df or df_train
        df_copy = df.copy().reset_index(drop=True)
        df_copy.loc[df_copy["Rxn"] , cols_to_f ] *= factor

        x_train = df_copy[features]
        y_train = df_copy[target]

        model = LinearRegression(fit_intercept=False)

        model.fit(x_train, y_train)
        if verbose : print("--------------------------------------------------------")
        if verbose : print("FACTOR =  ", factor)

        deltas_sum = 0.
        deltas_mean = 0.
        # calcolo le metriche sia di uno che dell'altro e vedo le differenze relative rispetto allo start
        for label, df_test in [ ("ABSOLUTE", df_abs), ("REACTION", df_rxn)]:
            # ("FULL",df),
            # df_test = df.copy().reset_index(drop=True) # To calculate metrics
            x_test = df_test[features]
            
            # Predizioni sul training
            y_pred = model.predict(x_test) + df_test["E_res_DZ"]
            y_true = df_test["E_QZ"]

            # Metriche
            rmsd_value = rmsd(y_true, y_pred)
            mape_value = mape_eps(y_true, y_pred, eps=1e-8)
            if label == ("ABSOLUTE"):
                deltar_rmsd = (rmsd_value - rmsd_value_abs)/rmsd_value_abs
                deltar_mape = (mape_value - mape_abs)/mape_abs
                plot_rmsd_abs.append(deltar_rmsd)
                plot_mape_abs.append(deltar_mape)
                if verbose : 
                    print("\nd_RMSD di ABSOLUTE = ", deltar_rmsd)
                    print("d_MAPE di ABSOLUTE = ", deltar_mape)
            if label == ("REACTION"):
                deltar_rmsd = (rmsd_value - rmsd_value_rxn)/rmsd_value_rxn
                deltar_mape = (mape_value - mape_rxn)/mape_rxn
                plot_rmsd_rxn.append(deltar_rmsd)
                plot_mape_rxn.append(deltar_mape)
                to_top_rmse_rxn = (rmsd_value - best_rmsd_value_rxn)/best_rmsd_value_rxn    
                if verbose : 
                    print("d_RMSD di REACTION = ", deltar_rmsd)
                    print("d_MAPE di REACTION = ", deltar_mape)
                    print ("to top : ", to_top_rmse_rxn)
            # if deltar_rmsd*deltar_mape > 0:
            if deltar_rmsd < 0 and deltar_mape < 0:
                deltas_sum = deltas_sum + deltar_rmsd + deltar_mape
            else :
                deltas_sum = np.inf
            
            
        
        deltas_mean = deltas_sum / 4 

        if verbose : print("deltas mean tra i quattro", deltas_mean)

        if deltas_mean < best_deltas and to_top_rmse_rxn < 0.05:
            best_deltas = deltas_mean
            best_factor = factor
            best_betas = {
                'coefs': dict(zip(features, model.coef_)),
                'intercept': model.intercept_
            }

        if verbose : print("best deltas", best_deltas)
        
    if best_factor == None:
        print("\n!!! ERROR !!!")
        print("no best_factor founded")
        break

    print("--------------------------------------------------------")
    # GIRO FINALE CON IL BEST FACTOR

    # best_factor = 23.604
    #print("\nBest Factor: ", best_factor)
    #print("\nBest Coefs: ", best_betas)

    df_copy = df.copy().reset_index(drop=True)
    df_copy.loc[df_copy["Rxn"] , cols_to_f] *= best_factor
    x_train = df_copy[features]
    y_train = df_copy[target]

    model = LinearRegression(fit_intercept=False)

    model.fit(x_train, y_train)


    for label, df_test in [ ("ABSOLUTE", df_abs), ("REACTION", df_rxn)]:

        x_test = df_test[features]
        y_pred = model.predict(x_test) + df_test["E_res_DZ"]
        y_true = df_test["E_QZ"]

        # Metriche
        rmsd_value = rmsd(y_true, y_pred)
        mape_value = mape_eps(y_true, y_pred, eps=1e-8)
        #print("\n  METRICS ",label, " ENERGY:")
        #print("RMSD :",rmsd_value)
        #print("MAPE :",mape_value)

    #print("Miglioramento finale delle metriche", best_deltas*100, "%")



#----------------------------------------------------
# PER FARE IL GRAFICO DEI DELTAS E SALVARE I PUNTI IN UN FILE plot.csv
#----------------------------------------------------

# Riordinare le varie colonne per poter ottenere un grafico senza sbaffo
# quando fa più di un ciclo
if g > 0:
    factors_arr = np.array(plot_factor)
    indices = np.argsort(factors_arr)
    plot_factor = factors_arr[indices]
    plot_rmsd_abs = np.array(plot_rmsd_abs)[indices]
    plot_rmsd_rxn = np.array(plot_rmsd_rxn)[indices]
    plot_mape_abs = np.array(plot_mape_abs)[indices]
    plot_mape_rxn = np.array(plot_mape_rxn)[indices]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# --- Grafico Delta RMSD ---
ax1.plot(plot_factor, plot_rmsd_abs, 'o-', label='Absolute', color='blue')
ax1.plot(plot_factor, plot_rmsd_rxn, 'o-', label='Reaction', color='red')
ax1.axhline(0, color='black', linestyle='--', alpha=0.5)
ax1.set_title("$\Delta$ RMSD Relativo %")
ax1.set_xlabel("Factor")
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.axvline(best_factor, color='black', linestyle='-', alpha=0.3)
# linee massimali calcolate in un check node
max_abs_rmsd = (best_rmsd_value_abs - rmsd_value_abs)/rmsd_value_abs
max_rxn_rmsd = (best_rmsd_value_rxn - rmsd_value_rxn)/rmsd_value_rxn
ax1.axhline(max_abs_rmsd, color='blue', linestyle='--', alpha=0.3)
ax1.axhline(max_rxn_rmsd, color='red', linestyle='--', alpha=0.3)

# --- Grafico Delta MAPE ---
ax2.plot(plot_factor, plot_mape_abs, 's-', label='Absolute', color='blue')
ax2.plot(plot_factor, plot_mape_rxn, 's-', label='Reaction', color='red')
ax2.axhline(0, color='black', linestyle='--', alpha=0.5)
ax2.set_title("$\Delta$ MAPE Relativo %")
ax2.set_xlabel("Factor")
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.axvline(best_factor, color='black', linestyle='-', alpha=0.3)
# linee massimali calcolate in un check node
max_abs_mape = (best_mape_value_abs - mape_abs)/mape_abs
max_rxn_mape = (best_mape_value_rxn - mape_rxn)/mape_rxn
ax2.axhline(max_abs_mape, color='blue', linestyle='--', alpha=0.3)
ax2.axhline(max_rxn_mape, color='red', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig("deltas.png")
# plt.show()



# SALVATAGGIO IN UN FILE .csv
if True:
    data_to_save = {
        'factor': plot_factor,
        'ABS_rmsd': plot_rmsd_abs,
        'ABS_mape': plot_mape_abs,
        'RXN_rmsd': plot_rmsd_rxn,
        'RXN_mape': plot_mape_rxn
    }
    df_panda = pd.DataFrame(data_to_save)
    df_panda.to_csv("deltas.csv", index=False, sep=';')


#----------------------------------------------------
# GIRO FINALE CON IL BEST FACTOR 
# O FACTOR A SCELTA
#----------------------------------------------------


# best_factor = 34.0
print("\nBest Factor: ", best_factor)

df_copy = df.copy().reset_index(drop=True)
df_copy.loc[df_copy["Rxn"] , cols_to_f] *= best_factor
x_train = df_copy[features]
y_train = df_copy[target]

model = LinearRegression(fit_intercept=False)

model.fit(x_train, y_train)

deltas_sum = 0.
for label, df_test in [ ("ABSOLUTE", df_abs), ("REACTION", df_rxn)]:

    x_test = df_test[features]
    y_pred = model.predict(x_test) + df_test["E_res_DZ"]
    y_true = df_test["E_QZ"]

    # Metriche
    rmsd_value = rmsd(y_true, y_pred)
    mape_value = mape_eps(y_true, y_pred, eps=1e-8)
    print("\n  METRICS ",label, " ENERGY:")
    print("RMSD :",rmsd_value)
    print("MAPE :",mape_value)

    if label == ("ABSOLUTE"):
        deltar_rmsd = (rmsd_value - rmsd_value_abs)/rmsd_value_abs
        deltar_mape = (mape_value - mape_abs)/mape_abs
        print("\nd_RMSD di ABSOLUTE = ", deltar_rmsd)
        print("d_MAPE di ABSOLUTE = ", deltar_mape)
    if label == ("REACTION"):
        deltar_rmsd = (rmsd_value - rmsd_value_rxn)/rmsd_value_rxn
        deltar_mape = (mape_value - mape_rxn)/mape_rxn
        print("\nd_RMSD di REACTION = ", deltar_rmsd)
        print("d_MAPE di REACTION = ", deltar_mape)

    deltas_sum = deltas_sum + deltar_rmsd + deltar_mape
    deltas_mean = deltas_sum / 4 


print("Miglioramento finale delle metriche", deltas_mean*100, "%")
print("\n", dict(zip(features, model.coef_)),
                'intercept:', model.intercept_)

print("\nFACTOR GRID SEARCH SAVED IN deltas.png & deltas.csv")