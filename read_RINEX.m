clear all
close all
clc

gnss_log = readtable('raw.csv');
eph = navsu.readfiles.loadRinexNav('BRDC00WRD_R_20210590000_01D_CN.rnx');
eph_EN = navsu.readfiles.loadRinexNav('BRDC00WRD_R_20210590000_01D_EN.rnx');
eph_GN = navsu.readfiles.loadRinexNav('BRDC00WRD_R_20210590000_01D_GN.rnx');
eph_RN = navsu.readfiles.loadRinexNav('BRDC00WRD_R_20210590000_01D_RN.rnx');

eph.gal = eph_EN.gal;
eph.gps = eph_GN.gps;
eph.glo = eph_RN.glo;

