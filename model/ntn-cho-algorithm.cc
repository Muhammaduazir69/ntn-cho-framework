/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2026
 *
 * SPDX-License-Identifier: GPL-2.0-only
 *
 * Author: Muhammad Uzair
 */

#include "ntn-cho-algorithm.h"

#include <ns3/enum.h>
#include <ns3/log.h>
#include <ns3/simulator.h>

#include <algorithm>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE("NtnChoAlgorithm");
NS_OBJECT_ENSURE_REGISTERED(NtnChoAlgorithm);

TypeId
NtnChoAlgorithm::GetTypeId()
{
    static TypeId tid =
        TypeId("ns3::NtnChoAlgorithm")
            .SetParent<Object>()
            .SetGroupName("NtnCho")
            .AddConstructor<NtnChoAlgorithm>()
            .AddTraceSource("HandoverExecuted",
                            "Fired when a CHO handover is executed",
                            MakeTraceSourceAccessor(&NtnChoAlgorithm::m_handoverExecutedTrace),
                            "ns3::NtnChoAlgorithm::HandoverExecutedTracedCallback")
            .AddTraceSource("HandoverOutcome",
                            "Fired with handover success/failure outcome",
                            MakeTraceSourceAccessor(&NtnChoAlgorithm::m_handoverOutcomeTrace),
                            "ns3::NtnChoAlgorithm::HandoverOutcomeTracedCallback")
            .AddTraceSource("CandidateEval",
                            "Fired when a candidate is evaluated",
                            MakeTraceSourceAccessor(&NtnChoAlgorithm::m_candidateEvalTrace),
                            "ns3::NtnChoAlgorithm::CandidateEvalTracedCallback")
            .AddTraceSource("StateTransition",
                            "Fired on CHO state machine transition",
                            MakeTraceSourceAccessor(&NtnChoAlgorithm::m_stateTransitionTrace),
                            "ns3::NtnChoAlgorithm::StateTransitionTracedCallback");
    return tid;
}

NtnChoAlgorithm::NtnChoAlgorithm()
    : m_state(CHO_IDLE),
      m_servingCellId(INVALID_CELL_ID),
      m_ueVelocity(Vector(0, 0, 0)),
      m_lastHoTime(Seconds(0)),
      m_lastSourceCell(INVALID_CELL_ID)
{
    NS_LOG_FUNCTION(this);
}

NtnChoAlgorithm::~NtnChoAlgorithm()
{
    NS_LOG_FUNCTION(this);
}

void
NtnChoAlgorithm::DoDispose()
{
    NS_LOG_FUNCTION(this);
    StopMonitoring();
    m_tteEstimator = nullptr;
    m_orbitPredictor = nullptr;
    m_candidates.clear();
    Object::DoDispose();
}

void
NtnChoAlgorithm::Configure(ChoConfig config)
{
    NS_LOG_FUNCTION(this);
    m_config = config;
    NS_LOG_INFO("CHO configured: trigger=" << config.triggerType
                << " d1=" << config.d1Threshold_m
                << "m qualityTh=" << config.qualityThreshold_dB
                << "dB tteMin=" << config.tteMinimum.GetSeconds() << "s");
}

NtnChoAlgorithm::ChoConfig
NtnChoAlgorithm::GetConfig() const
{
    return m_config;
}

void
NtnChoAlgorithm::SetTteEstimator(Ptr<NtnTteEstimator> estimator)
{
    m_tteEstimator = estimator;
}

void
NtnChoAlgorithm::SetOrbitPredictor(Ptr<NtnOrbitPredictor> predictor)
{
    m_orbitPredictor = predictor;
}

void
NtnChoAlgorithm::AddCandidateCell(uint16_t cellId, uint32_t satId, uint32_t beamId)
{
    NS_LOG_FUNCTION(this << cellId << satId << beamId);

    if (m_candidates.size() >= m_config.maxCandidates)
    {
        NS_LOG_WARN("Maximum candidates (" << (int)m_config.maxCandidates
                    << ") reached. Ignoring cell " << cellId);
        return;
    }

    CandidateInfo cand;
    cand.cellId = cellId;
    cand.satId = satId;
    cand.beamId = beamId;
    m_candidates[cellId] = cand;

    if (m_state == CHO_IDLE)
    {
        TransitionState(CHO_PREPARED);
    }
}

void
NtnChoAlgorithm::RemoveCandidateCell(uint16_t cellId)
{
    m_candidates.erase(cellId);
}

void
NtnChoAlgorithm::ClearCandidates()
{
    m_candidates.clear();
    if (m_state == CHO_PREPARED || m_state == CHO_CONDITION_MONITORING)
    {
        StopMonitoring();
        TransitionState(CHO_IDLE);
    }
}

void
NtnChoAlgorithm::UpdateMeasurement(uint16_t cellId, double sinr_dB, double gain_dB)
{
    auto it = m_candidates.find(cellId);
    if (it != m_candidates.end())
    {
        it->second.sinr_dB = sinr_dB;
        it->second.gain_dB = gain_dB;
        it->second.lastUpdate = Simulator::Now();
    }
}

void
NtnChoAlgorithm::StartMonitoring(GeoCoordinate uePosition, Vector ueVelocity)
{
    NS_LOG_FUNCTION(this);
    m_uePosition = uePosition;
    m_ueVelocity = ueVelocity;

    if (m_state == CHO_PREPARED || m_state == CHO_CONDITION_MONITORING)
    {
        TransitionState(CHO_CONDITION_MONITORING);
        DoEvaluateConditions();
    }
}

void
NtnChoAlgorithm::StopMonitoring()
{
    NS_LOG_FUNCTION(this);
    if (m_monitorEvent.IsPending())
    {
        Simulator::Cancel(m_monitorEvent);
    }
}

void
NtnChoAlgorithm::DoEvaluateConditions()
{
    if (m_state != CHO_CONDITION_MONITORING)
    {
        return;
    }

    EvaluateConditions();

    // Schedule next evaluation
    m_monitorEvent = Simulator::Schedule(m_config.conditionMonitorPeriod,
                                         &NtnChoAlgorithm::DoEvaluateConditions,
                                         this);
}

void
NtnChoAlgorithm::EvaluateConditions()
{
    NS_LOG_FUNCTION(this);

    for (auto& [cellId, cand] : m_candidates)
    {
        // Step 1: Check D1 condition
        bool d1Now = CheckD1Condition(cand);

        if (d1Now && !cand.d1Met)
        {
            cand.d1Met = true;
            cand.d1MetSince = Simulator::Now();
            NS_LOG_DEBUG("D1 condition MET for cell " << cellId);
        }
        else if (!d1Now)
        {
            cand.d1Met = false;
            cand.admitted = false;
        }

        // Step 2: If D1 met and quality sufficient, compute TTE
        if (cand.d1Met && cand.sinr_dB >= m_config.qualityThreshold_dB)
        {
            if (m_tteEstimator && m_config.triggerType == TRIGGER_TTE_AWARE)
            {
                auto tteResult = m_tteEstimator->ComputeTte(
                    m_uePosition, m_ueVelocity,
                    cand.satId, cand.beamId,
                    m_config.gainThreshold_dB);

                cand.tte = tteResult.tte;

                // Admit if TTE >= minimum
                cand.admitted = (cand.tte >= m_config.tteMinimum);
            }
            else if (m_config.triggerType == TRIGGER_LOCATION_D1)
            {
                // Baseline: admit based on D1 + quality only (no TTE)
                cand.admitted = true;
                cand.tte = Seconds(0); // Unknown
            }
            else if (m_config.triggerType == TRIGGER_EVENT_A3)
            {
                // A3: admit based purely on SINR comparison
                cand.admitted = true;
                cand.tte = Seconds(0);
            }
        }
        else if (cand.d1Met && cand.sinr_dB < m_config.qualityThreshold_dB)
        {
            cand.admitted = false;
        }

        // Fire trace
        m_candidateEvalTrace(cellId, cand.sinr_dB, cand.tte, cand.admitted);

        if (cand.admitted && !m_admitCallback.IsNull())
        {
            m_admitCallback(cellId, cand.sinr_dB, cand.tte);
        }
    }
}

uint16_t
NtnChoAlgorithm::SelectBestCandidate() const
{
    NS_LOG_FUNCTION(this);

    // ================================================================
    // THE NOVEL TTE-AWARE CANDIDATE SELECTION ALGORITHM
    //
    // Input:  Set of candidate cells with measured SINR and computed TTE
    // Output: Best candidate cell ID
    //
    // 1. Filter: remove candidates with SINR < Q_threshold
    // 2. Filter: remove candidates with TTE < TTE_minimum
    // 3. Select: candidate with maximum TTE
    // 4. Tie-break: if multiple within epsilon, pick highest SINR
    // ================================================================

    std::vector<const CandidateInfo*> admissible;

    for (const auto& [cellId, info] : m_candidates)
    {
        if (info.admitted &&
            info.d1Met &&
            info.sinr_dB >= m_config.qualityThreshold_dB &&
            info.tte >= m_config.tteMinimum)
        {
            admissible.push_back(&info);
        }
    }

    if (admissible.empty())
    {
        NS_LOG_DEBUG("No admissible candidates for TTE-aware selection");
        return INVALID_CELL_ID;
    }

    // Sort by TTE descending
    std::sort(admissible.begin(), admissible.end(),
              [](const CandidateInfo* a, const CandidateInfo* b) {
                  return a->tte > b->tte;
              });

    // Tie-breaking: within epsilon of the best TTE, pick highest SINR
    const CandidateInfo* best = admissible[0];
    for (size_t i = 1; i < admissible.size(); i++)
    {
        if ((best->tte - admissible[i]->tte) <= m_config.tteEpsilon &&
            admissible[i]->sinr_dB > best->sinr_dB)
        {
            best = admissible[i];
        }
    }

    NS_LOG_INFO("TTE-aware selection: cell " << best->cellId
                << " (TTE=" << best->tte.GetSeconds()
                << "s, SINR=" << best->sinr_dB << " dB)");
    return best->cellId;
}

uint16_t
NtnChoAlgorithm::SelectBaselineA3(double servingSinr_dB) const
{
    NS_LOG_FUNCTION(this << servingSinr_dB);

    uint16_t bestCell = INVALID_CELL_ID;
    double bestSinr = -200.0;

    for (const auto& [cellId, info] : m_candidates)
    {
        // A3: neighbor SINR > serving SINR + offset
        if (info.sinr_dB > (servingSinr_dB + m_config.a3Offset_dB) &&
            info.sinr_dB > bestSinr)
        {
            bestSinr = info.sinr_dB;
            bestCell = cellId;
        }
    }

    return bestCell;
}

uint16_t
NtnChoAlgorithm::SelectBaselineLocationOnly() const
{
    NS_LOG_FUNCTION(this);

    uint16_t bestCell = INVALID_CELL_ID;
    double bestSinr = -200.0;

    for (const auto& [cellId, info] : m_candidates)
    {
        if (info.d1Met && info.sinr_dB >= m_config.qualityThreshold_dB &&
            info.sinr_dB > bestSinr)
        {
            bestSinr = info.sinr_dB;
            bestCell = cellId;
        }
    }

    return bestCell;
}

void
NtnChoAlgorithm::ExecuteHandover(uint16_t targetCellId)
{
    NS_LOG_FUNCTION(this << targetCellId);

    if (targetCellId == INVALID_CELL_ID)
    {
        NS_LOG_WARN("Cannot execute HO to invalid cell");
        m_handoverOutcomeTrace(targetCellId, false, "InvalidTarget");
        return;
    }

    // Record Time-of-Stay from last handover
    Time tos = Seconds(0);
    if (m_lastHoTime > Seconds(0))
    {
        tos = Simulator::Now() - m_lastHoTime;
    }

    uint16_t sourceCell = m_servingCellId;

    // Check for ping-pong
    bool isPingPong = (targetCellId == m_lastSourceCell &&
                       tos < Seconds(5.0)); // 5s ping-pong threshold

    TransitionState(CHO_EXECUTING);

    // Fire execution trace
    m_handoverExecutedTrace(sourceCell, targetCellId, tos);

    // Start T304 timer
    m_t304Event = Simulator::Schedule(m_config.t304Timer, [this, targetCellId]() {
        // T304 expired - handover failure
        NS_LOG_WARN("T304 expired for HO to cell " << targetCellId);
        m_handoverOutcomeTrace(targetCellId, false, "T304Expired");
        TransitionState(CHO_IDLE);
    });

    // Execute via callback
    if (!m_hoCallback.IsNull())
    {
        m_hoCallback(sourceCell, targetCellId);
    }

    // Update state
    m_lastSourceCell = sourceCell;
    m_lastHoTime = Simulator::Now();
    m_servingCellId = targetCellId;

    // Cancel T304 (successful execution)
    if (m_t304Event.IsPending())
    {
        Simulator::Cancel(m_t304Event);
    }

    if (isPingPong)
    {
        m_handoverOutcomeTrace(targetCellId, true, "PingPong");
    }
    else
    {
        m_handoverOutcomeTrace(targetCellId, true, "Success");
    }

    TransitionState(CHO_COMPLETED);

    // Reset to monitoring for next CHO cycle
    m_candidates.clear();
    TransitionState(CHO_IDLE);
}

void
NtnChoAlgorithm::CancelHandover()
{
    NS_LOG_FUNCTION(this);
    StopMonitoring();
    if (m_t304Event.IsPending())
    {
        Simulator::Cancel(m_t304Event);
    }
    TransitionState(CHO_IDLE);
}

NtnChoAlgorithm::ChoState
NtnChoAlgorithm::GetState() const
{
    return m_state;
}

std::map<uint16_t, NtnChoAlgorithm::CandidateInfo>
NtnChoAlgorithm::GetCandidates() const
{
    return m_candidates;
}

uint32_t
NtnChoAlgorithm::GetNumAdmittedCandidates() const
{
    uint32_t count = 0;
    for (const auto& [id, info] : m_candidates)
    {
        if (info.admitted)
        {
            count++;
        }
    }
    return count;
}

void
NtnChoAlgorithm::SetHandoverExecutionCallback(HandoverExecutionCallback cb)
{
    m_hoCallback = cb;
}

void
NtnChoAlgorithm::SetCandidateAdmittedCallback(CandidateAdmittedCallback cb)
{
    m_admitCallback = cb;
}

bool
NtnChoAlgorithm::CheckD1Condition(const CandidateInfo& cand) const
{
    if (!m_orbitPredictor)
    {
        return false;
    }

    auto snap = m_orbitPredictor->GetBeamSnapshot(cand.satId, cand.beamId, m_uePosition);

    // D1: UE within d1Threshold of beam center
    Vector ueCart = m_uePosition.ToVector();
    Vector bcCart = snap.beamCenter.ToVector();
    double dx = ueCart.x - bcCart.x;
    double dy = ueCart.y - bcCart.y;
    double dz = ueCart.z - bcCart.z;
    double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    return (dist <= m_config.d1Threshold_m);
}

void
NtnChoAlgorithm::TransitionState(ChoState newState)
{
    if (m_state != newState)
    {
        NS_LOG_INFO("CHO state: " << m_state << " -> " << newState);
        ChoState old = m_state;
        m_state = newState;
        m_stateTransitionTrace(old, newState);
    }
}

} // namespace ns3
