# %% Import modules to work
import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.preprocessing import PolynomialFeatures

# %% Creation of dictionary with gene length
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cya_clusterPipeline')
WH7803_genome = SeqIO.parse('/Users/u1866168/Documents/OneDrive - University of '
                            'Warwick/Genome_Database/Cya_clusterPipeline/23615-Syn_WH7803.gb', 'genbank')
Gen_len_Dict = {}

for genome in WH7803_genome:
    for feature in genome.features:
        if feature.type == "CDS":  # Apply to those that are not manually created and might have CyaClusterID
            ID = str(feature.qualifiers["locus_tag"][0])
            length = len(str(feature.extract(genome.seq)))
            Gen_len_Dict.update({ID: length})
# %% Preparation of Temporal models Plus P Non infected
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Branko_Host_analysis')
Dataframe = {}
Group1_ArrayX1 = np.array([0, 3, 6, 9])

Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P+_NoInf', index_col='Locustag')

for index, row in Data1.iterrows():
    #  Estimation of rpkm per gene and creation of Y array.
    Serie1 = [
        (int(row['0h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['3h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['6h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['9h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000)))]
    Serie2 = [
        (int(row['0h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['3h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['6h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['9h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000)))]
    Serie3 = [
        (int(row['0h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['3h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['6h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['9h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000)))]

    Group1_Arrayx1 = Group1_ArrayX1[:, np.newaxis]  # Reshape of X axis

    polynomial_features = PolynomialFeatures(degree=2)  # Extraction of features for linear model
    x_poly = polynomial_features.fit_transform(Group1_Arrayx1)  # Transformation of X axis
    model1 = LinearRegression()  # Setting the linear model class object
    model2 = LinearRegression()
    model3 = LinearRegression()
    model1.fit(x_poly, Serie1)  # Creation of the fit model for Data array. (2nd degree)
    model2.fit(x_poly, Serie2)
    model3.fit(x_poly, Serie3)
    y1_poly_pred = model1.predict(x_poly)  # Creation of Y data based on predicted model for R score
    y2_poly_pred = model2.predict(x_poly)
    y3_poly_pred = model3.predict(x_poly)

    # Extraction of coefficient, intercept and r score
    Dataframe.update({index: [model1.coef_[0], model1.coef_[1], model1.coef_[2], model1.intercept_,
                              r2_score(Serie1, y1_poly_pred),
                              model2.coef_[0], model2.coef_[1], model2.coef_[2], model2.intercept_,
                              r2_score(Serie2, y2_poly_pred),
                              model3.coef_[0], model3.coef_[1], model3.coef_[2], model3.intercept_,
                              r2_score(Serie3, y3_poly_pred)]})

# Preparation of Temporal model Minus P Non infected

Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P-_NoInf', index_col='Locustag')

for index, row in Data1.iterrows():
    #  Estimation of rpkm per gene and creation of Y array.
    Serie1 = [
        (int(row['0h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['3h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['6h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['9h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000)))]
    Serie2 = [
        (int(row['0h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['3h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['6h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['9h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000)))]
    Serie3 = [
        (int(row['0h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['3h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['6h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['9h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000)))]

    Group1_Arrayx1 = Group1_ArrayX1[:, np.newaxis]  # Reshape of X axis

    polynomial_features = PolynomialFeatures(degree=2)  # Extraction of features for linear model
    x_poly = polynomial_features.fit_transform(Group1_Arrayx1)  # Transformation of X axis
    model1 = LinearRegression()  # Setting the linear model class object
    model2 = LinearRegression()
    model3 = LinearRegression()
    model1.fit(x_poly, Serie1)  # Creation of the fit model for Data array. (2nd degree)
    model2.fit(x_poly, Serie2)
    model3.fit(x_poly, Serie3)
    y1_poly_pred = model1.predict(x_poly)  # Creation of Y data based on predicted model for R score
    y2_poly_pred = model2.predict(x_poly)
    y3_poly_pred = model3.predict(x_poly)

    Dataframe[index].append(model1.coef_[0])
    Dataframe[index].append(model1.coef_[1])
    Dataframe[index].append(model1.coef_[2])
    Dataframe[index].append(model1.intercept_)
    Dataframe[index].append(r2_score(Serie1, y1_poly_pred))
    Dataframe[index].append(model2.coef_[0])
    Dataframe[index].append(model2.coef_[1])
    Dataframe[index].append(model2.coef_[2])
    Dataframe[index].append(model2.intercept_)
    Dataframe[index].append(r2_score(Serie2, y2_poly_pred))
    Dataframe[index].append(model3.coef_[0])
    Dataframe[index].append(model3.coef_[1])
    Dataframe[index].append(model3.coef_[2])
    Dataframe[index].append(model3.intercept_)
    Dataframe[index].append(r2_score(Serie3, y3_poly_pred))

# Preparation of the temporal model for Plus P Infected

Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P+_Inf', index_col='Locustag')

for index, row in Data1.iterrows():
    #  Estimation of rpkm per gene and creation of Y array.
    Serie1 = [
        (int(row['0h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['3h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['6h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['9h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000)))]
    Serie2 = [
        (int(row['0h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['3h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['6h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['9h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000)))]
    Serie3 = [
        (int(row['0h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['3h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['6h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['9h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000)))]

    Group1_Arrayx1 = Group1_ArrayX1[:, np.newaxis]  # Reshape of X axis

    polynomial_features = PolynomialFeatures(degree=2)  # Extraction of features for linear model
    x_poly = polynomial_features.fit_transform(Group1_Arrayx1)  # Transformation of X axis
    model1 = LinearRegression()  # Setting the linear model class object
    model2 = LinearRegression()
    model3 = LinearRegression()
    model1.fit(x_poly, Serie1)  # Creation of the fit model for Data array. (2nd degree)
    model2.fit(x_poly, Serie2)
    model3.fit(x_poly, Serie3)
    y1_poly_pred = model1.predict(x_poly)  # Creation of Y data based on predicted model for R score
    y2_poly_pred = model2.predict(x_poly)
    y3_poly_pred = model3.predict(x_poly)

    Dataframe[index].append(model1.coef_[0])
    Dataframe[index].append(model1.coef_[1])
    Dataframe[index].append(model1.coef_[2])
    Dataframe[index].append(model1.intercept_)
    Dataframe[index].append(r2_score(Serie1, y1_poly_pred))
    Dataframe[index].append(model2.coef_[0])
    Dataframe[index].append(model2.coef_[1])
    Dataframe[index].append(model2.coef_[2])
    Dataframe[index].append(model2.intercept_)
    Dataframe[index].append(r2_score(Serie2, y2_poly_pred))
    Dataframe[index].append(model3.coef_[0])
    Dataframe[index].append(model3.coef_[1])
    Dataframe[index].append(model3.coef_[2])
    Dataframe[index].append(model3.intercept_)
    Dataframe[index].append(r2_score(Serie3, y3_poly_pred))

# Preparation of the temporal model for Minus P Infected

Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P-_Inf', index_col='Locustag')

for index, row in Data1.iterrows():
    #  Estimation of rpkm per gene and creation of Y array.
    Serie1 = [
        (int(row['0h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),

        (int(row['6h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000))),
        (int(row['9h_1']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_1'].sum() / 1000000))) / (int(row['0h_1']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_1'].sum() / 1000000)))]
    Serie2 = [
        (int(row['0h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['3h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['6h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000))),
        (int(row['9h_2']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_2'].sum() / 1000000))) / (int(row['0h_2']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_2'].sum() / 1000000)))]
    Serie3 = [
        (int(row['0h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['0h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['3h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['3h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['6h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['6h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000))),
        (int(row['9h_3']) / ((Gen_len_Dict[index] / 1000) * (Data1['9h_3'].sum() / 1000000))) / (int(row['0h_3']) /
                                                                                                 ((Gen_len_Dict[
                                                                                                       index] / 1000) * (
                                                                                                          Data1[
                                                                                                              '0h_3'].sum() / 1000000)))]

    Group1_Arrayx1 = Group1_ArrayX1[:, np.newaxis]
    Group1_Array_Serie1 = np.array([0, 6, 9])
    Group1_Array_Serie1 = Group1_Array_Serie1[:, np.newaxis]

    polynomial_features = PolynomialFeatures(degree=2)
    x_poly = polynomial_features.fit_transform(Group1_Arrayx1)
    xSerie1 = polynomial_features.fit_transform(Group1_Array_Serie1)
    model1 = LinearRegression()
    model2 = LinearRegression()
    model3 = LinearRegression()
    model1.fit(xSerie1, Serie1)
    model2.fit(x_poly, Serie2)
    model3.fit(x_poly, Serie3)
    y1_poly_pred = model1.predict(xSerie1)
    y2_poly_pred = model2.predict(x_poly)
    y3_poly_pred = model3.predict(x_poly)

    Dataframe[index].append(model1.coef_[0])
    Dataframe[index].append(model1.coef_[1])
    Dataframe[index].append(model1.coef_[2])
    Dataframe[index].append(model1.intercept_)
    Dataframe[index].append(r2_score(Serie1, y1_poly_pred))
    Dataframe[index].append(model2.coef_[0])
    Dataframe[index].append(model2.coef_[1])
    Dataframe[index].append(model2.coef_[2])
    Dataframe[index].append(model2.intercept_)
    Dataframe[index].append(r2_score(Serie2, y2_poly_pred))
    Dataframe[index].append(model3.coef_[0])
    Dataframe[index].append(model3.coef_[1])
    Dataframe[index].append(model3.coef_[2])
    Dataframe[index].append(model3.intercept_)
    Dataframe[index].append(r2_score(Serie3, y3_poly_pred))

TempDataframe = pd.DataFrame.from_dict(Dataframe, orient='index_genome')
# TempDataframe.to_csv('Temporal_model_Host.csv')
