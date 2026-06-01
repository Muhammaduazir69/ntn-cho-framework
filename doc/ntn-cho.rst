..
   SPDX-License-Identifier: GPL-2.0-only
   Copyright (c) 2026 Muhammad Uzair and contributors

ntn-cho Module
==============

.. include:: replace.txt
.. highlight:: cpp

Overview
--------

The ``ntn-cho`` module provides a 3GPP Release-17 Conditional Handover
(CHO) framework for Non-Terrestrial Networks (NTN) in ns-3.  It adds an
orbit-informed **Time-to-Exit (TTE)** trigger to the standard Rel-17
CHO state machine and ships the classical A3 / location / time-based
selectors as baselines for comparison.

The module **consumes** the SGP4 orbit propagation and the Walker
constellation geometry provided by the SNS3 ``satellite`` module
(``SatSGP4MobilityModel`` + antenna-gain patterns); it does not
re-implement orbit propagation or a constellation generator of its own.

Model description
-----------------

The source code lives in ``contrib/ntn-cho/model`` and
``contrib/ntn-cho/helper``.

Design
~~~~~~

The module exposes the following public classes:

* ``NtnTteEstimator`` — Time-to-Exit estimator.  Propagates the
  candidate satellite forward with a coarse step, then refines the
  beam-exit instant with a **binary search**; the per-candidate result
  is cached.
* ``NtnChoAlgorithm`` — Rel-17 CHO state machine
  (``IDLE → PREPARE → MONITOR → EXEC``).  The ``TriggerType`` enum
  selects between ``TRIGGER_EVENT_A3``, ``TRIGGER_LOCATION_D1``
  (condEventD1), ``TRIGGER_TIME_BASED``, the proposed
  ``TRIGGER_TTE_AWARE``, and ``TRIGGER_THZ_BEAM_QUALITY``.
* ``NtnOrbitPredictor`` — wraps the ``satellite`` module's
  ``SatSGP4MobilityModel`` and antenna-gain-pattern container to
  predict satellite sub-point, elevation, and beam coverage over time.
* ``NtnMeasurementModel`` — per-UE RSRP / SINR / elevation / Doppler
  generation that feeds the CHO triggers.
* ``NtnAiInterface`` — optional ns3-ai bridge (observation / action
  structs) for learning-based HO policies.  Requires the ``ns3-ai``
  module.

Helpers:

* ``NtnChoHelper`` — wires the measurement model, TTE estimator, and
  CHO algorithm onto a node container.
* ``NtnRealisticMobilityHelper`` — per-class ground-UE motion for the
  seven TR 38.811 §6.1.1.1 UE classes.
* ``NtnRealisticTrafficHelper`` — per-class UDP traffic profiles.

Scope and limitations
~~~~~~~~~~~~~~~~~~~~~~~

* The module targets LEO NTN scenarios; GEO/MEO geometries are not
  validated.
* The CHO procedure follows TS 38.331 Release 17 conditional handover.
  Multi-connectivity (split-CHO) is not implemented.
* Orbit propagation and constellation geometry are delegated to the
  ``satellite`` module; ``ntn-cho`` does not own that code.

Attributes
~~~~~~~~~~

The registered, commonly-tuned attributes are:

* ``NtnTteEstimator::PredictionStep`` — coarse forward step.
* ``NtnTteEstimator::MaxPredictionWindow`` — forward search horizon.
* ``NtnTteEstimator::BinarySearchTolerance`` — exit-time precision
  (default 100 ms).
* ``NtnTteEstimator::MaxBinarySearchIterations``.
* ``NtnOrbitPredictor::MinGainThreshold`` — beam-edge gain threshold.
* ``NtnChoAlgorithm::EnableMultiBandCho``,
  ``ThzBeamTrackingThreshold``, ``ThzSnrThreshold``,
  ``ThzBeamwidth`` — multi-band / THz-beam options.

The remaining trigger parameters (D1 distance threshold, A3 offset and
time-to-trigger, SINR quality threshold, minimum TTE) are set per run
via the example command-line arguments
(``--d1Threshold``, ``--qualityTh``, ``--tteMinimum``, …).

Output
~~~~~~

``ntn-cho-full-constellation`` writes to ``--outputDir``:
``handover_events.csv``, ``measurements.csv``, ``tte_computations.csv``,
``satellite_tracks.csv``, ``ue_tracks.csv``, ``kpi_timeseries.csv``,
the matching GeoJSON layers (``satellite_positions.geojson``,
``ue_positions.geojson``, ``beam_footprints.geojson``,
``handover_events.geojson``), and ``kpi_summary.txt``.

Examples
~~~~~~~~

* ``ntn-cho-leo-basic.cc`` — minimal single-satellite CHO walkthrough.
* ``ntn-cho-full-constellation.cc`` — multi-satellite Walker
  constellation; selects the algorithm with
  ``--algorithm=a3|location|time|tte-aware``.
* ``ntn-realistic-mobility-demo.cc`` — seven-class UE mobility demo.

Testing
-------

The ``ntn-cho`` unit-test suite (``test/ntn-cho-test-suite.cc``)
contains three ``QUICK`` cases — ``NtnChoAlgorithmTestCase``,
``NtnChoStateMachineTestCase``, and ``NtnMeasurementModelTestCase`` —
which validate the algorithm and state-machine logic:

.. sourcecode:: bash

   ./test.py --suite=ntn-cho

The aggregate KPIs reported in the associated manuscript (handover
counts, ping-pong rate, success rate) are produced by sweeping
``ntn-cho-full-constellation`` over multiple ``--rngRun`` seeds and the
four ``--algorithm`` settings; they are **not** asserted by the unit
suite.

References
~~~~~~~~~~

* 3GPP TS 38.331 v17, Radio Resource Control (RRC); Protocol
  specification.
* 3GPP TR 38.821 v16, Solutions for NR to support NTN.
* 3GPP TR 38.811 v15, Study on NR to support NTN.
