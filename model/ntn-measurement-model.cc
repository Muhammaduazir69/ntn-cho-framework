/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2026
 *
 * SPDX-License-Identifier: GPL-2.0-only
 *
 * Author: Muhammad Uzair
 */

#include "ntn-measurement-model.h"

#include <ns3/log.h>
#include <ns3/satellite-mobility-model.h>
#include <ns3/simulator.h>

#include <algorithm>
#include <cmath>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE("NtnMeasurementModel");
NS_OBJECT_ENSURE_REGISTERED(NtnMeasurementModel);

TypeId
NtnMeasurementModel::GetTypeId()
{
    static TypeId tid =
        TypeId("ns3::NtnMeasurementModel")
            .SetParent<Object>()
            .SetGroupName("NtnCho")
            .AddConstructor<NtnMeasurementModel>()
            .AddTraceSource("MeasurementReport",
                            "Fired when a measurement is computed",
                            MakeTraceSourceAccessor(&NtnMeasurementModel::m_measurementTrace),
                            "ns3::NtnMeasurementModel::MeasurementTracedCallback");
    return tid;
}

NtnMeasurementModel::NtnMeasurementModel()
    : m_scenario(NTN_SUBURBAN),
      m_carrierFreqHz(2.0e9),   // S-band default
      m_bandwidthHz(30.0e6),    // 30 MHz default
      m_satTxPower_dBm(40.0),   // 10W per beam
      m_ueNoiseFigure_dB(7.0)
{
    NS_LOG_FUNCTION(this);
}

NtnMeasurementModel::~NtnMeasurementModel()
{
    NS_LOG_FUNCTION(this);
}

void
NtnMeasurementModel::DoDispose()
{
    m_orbitPredictor = nullptr;
    Object::DoDispose();
}

void
NtnMeasurementModel::SetOrbitPredictor(Ptr<NtnOrbitPredictor> predictor)
{
    m_orbitPredictor = predictor;
}

void
NtnMeasurementModel::SetNtnScenario(NtnScenario scenario)
{
    m_scenario = scenario;
}

void
NtnMeasurementModel::SetCarrierFrequency(double freqHz)
{
    m_carrierFreqHz = freqHz;
}

void
NtnMeasurementModel::SetBandwidth(double bwHz)
{
    m_bandwidthHz = bwHz;
}

void
NtnMeasurementModel::SetSatelliteTxPower(double txPower_dBm)
{
    m_satTxPower_dBm = txPower_dBm;
}

void
NtnMeasurementModel::SetUeNoiseFigure(double nf_dB)
{
    m_ueNoiseFigure_dB = nf_dB;
}

NtnMeasurementModel::Measurement
NtnMeasurementModel::ComputeMeasurement(GeoCoordinate uePosition,
                                         uint32_t satId,
                                         uint32_t beamId) const
{
    NS_LOG_FUNCTION(this << satId << beamId);
    NS_ASSERT_MSG(m_orbitPredictor, "OrbitPredictor not set");

    Measurement meas;
    meas.satId = satId;
    meas.beamId = beamId;

    // Get beam snapshot
    auto snap = m_orbitPredictor->GetBeamSnapshot(satId, beamId, uePosition);
    meas.elevationAngle_deg = snap.elevationAngle_deg;
    meas.slantRange_km = snap.slantRange_km;
    meas.propagationDelay = snap.propagationDelay;
    meas.antennaGain_dB = snap.gainAtUe_dB;

    // Compute path loss: Free-space + atmospheric/clutter (3GPP TR 38.811)
    double distance_m = snap.slantRange_km * 1000.0;
    double fspl_dB = ComputeFreeSpacePathLoss(distance_m, m_carrierFreqHz);

    // Additional NTN-specific losses per TR 38.811
    // Atmospheric absorption (simplified): ~0.07 dB/km for S-band at 10deg elevation
    double atmosphericLoss_dB = 0.0;
    if (meas.elevationAngle_deg < 20.0)
    {
        atmosphericLoss_dB = 0.5; // Increased loss at low elevation
    }

    // Scintillation loss (simplified): ~0.5 dB for S-band
    double scintillationLoss_dB = 0.5;

    // Clutter loss depends on scenario (TR 38.811 Table 6.6.2-1)
    double clutterLoss_dB = 0.0;
    switch (m_scenario)
    {
    case NTN_DENSE_URBAN:
        clutterLoss_dB = 4.0;
        break;
    case NTN_URBAN:
        clutterLoss_dB = 2.5;
        break;
    case NTN_SUBURBAN:
        clutterLoss_dB = 1.0;
        break;
    case NTN_RURAL:
        clutterLoss_dB = 0.5;
        break;
    }

    // Shadow fading margin (depends on LoS probability which depends on elevation)
    double shadowFading_dB = 0.0;
    if (meas.elevationAngle_deg < 30.0)
    {
        shadowFading_dB = 3.0; // Larger fading at low elevation
    }
    else if (meas.elevationAngle_deg < 60.0)
    {
        shadowFading_dB = 2.0;
    }
    else
    {
        shadowFading_dB = 1.0;
    }

    meas.pathLoss_dB = fspl_dB + atmosphericLoss_dB + scintillationLoss_dB +
                       clutterLoss_dB + shadowFading_dB;

    // RSRP = TxPower + AntennaGain - PathLoss
    meas.rsrp_dBm = m_satTxPower_dBm + meas.antennaGain_dB - meas.pathLoss_dB;

    // SINR = RSRP - NoisePower (simplified: single-link, no interference)
    double noisePower_dBm = ComputeNoisePower();
    meas.sinr_dB = meas.rsrp_dBm - noisePower_dBm;

    // Compute Doppler shift
    meas.dopplerShift_Hz = ComputeDopplerShift(uePosition, satId);

    NS_LOG_DEBUG("Measurement: sat=" << satId << " beam=" << beamId
                 << " RSRP=" << meas.rsrp_dBm << " dBm"
                 << " SINR=" << meas.sinr_dB << " dB"
                 << " elev=" << meas.elevationAngle_deg << " deg"
                 << " range=" << meas.slantRange_km << " km");

    m_measurementTrace(satId, beamId, meas.rsrp_dBm, meas.sinr_dB);
    return meas;
}

std::vector<NtnMeasurementModel::Measurement>
NtnMeasurementModel::ScanVisibleBeams(GeoCoordinate uePosition,
                                       double minElevation_deg) const
{
    NS_LOG_FUNCTION(this << minElevation_deg);

    std::vector<Measurement> results;

    auto visible = m_orbitPredictor->GetVisibleSatellites(uePosition, minElevation_deg);

    uint16_t cellIdCounter = 1;
    for (const auto& vs : visible)
    {
        Measurement meas = ComputeMeasurement(uePosition, vs.satId, vs.bestBeamId);
        meas.cellId = cellIdCounter++;
        results.push_back(meas);
    }

    // Sort by SINR descending
    std::sort(results.begin(), results.end(),
              [](const Measurement& a, const Measurement& b) {
                  return a.sinr_dB > b.sinr_dB;
              });

    return results;
}

double
NtnMeasurementModel::ComputeDopplerShift(GeoCoordinate uePosition, uint32_t satId) const
{
    NS_LOG_FUNCTION(this << satId);

    Ptr<SatMobilityModel> satMob = m_orbitPredictor->GetSatelliteMobility(satId);
    GeoCoordinate satPos = satMob->GetGeoPosition();

    // Get satellite velocity from mobility model
    Vector satVel = satMob->GetVelocity(); // m/s in ECEF

    // Vector from UE to satellite
    Vector ueCart = uePosition.ToVector();
    Vector satCart = satPos.ToVector();
    double dx = satCart.x - ueCart.x;
    double dy = satCart.y - ueCart.y;
    double dz = satCart.z - ueCart.z;
    double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (dist < 1.0)
    {
        return 0.0;
    }

    // Unit vector from UE to satellite
    double ux = dx / dist;
    double uy = dy / dist;
    double uz = dz / dist;

    // Radial velocity component (positive = approaching)
    double vRadial = -(satVel.x * ux + satVel.y * uy + satVel.z * uz);

    // Doppler shift: f_d = (v_radial / c) * f_c
    double dopplerHz = (vRadial / SPEED_OF_LIGHT) * m_carrierFreqHz;

    return dopplerHz;
}

double
NtnMeasurementModel::ComputeFreeSpacePathLoss(double distance_m, double freqHz) const
{
    // FSPL = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    double fspl = 20.0 * std::log10(distance_m) +
                  20.0 * std::log10(freqHz) +
                  20.0 * std::log10(4.0 * M_PI / SPEED_OF_LIGHT);
    return fspl;
}

double
NtnMeasurementModel::ComputeNoisePower() const
{
    // N = 10*log10(k*T*B) + NF
    // k*T at 290K = -174 dBm/Hz
    double thermalNoise_dBmHz = -174.0;
    double noisePower_dBm = thermalNoise_dBmHz + 10.0 * std::log10(m_bandwidthHz) +
                            m_ueNoiseFigure_dB;
    return noisePower_dBm;
}

} // namespace ns3
