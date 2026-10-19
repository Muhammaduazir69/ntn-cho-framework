build_lib_example(
  NAME ntn-cho-leo-basic
  SOURCE_FILES ntn-cho-leo-basic.cc
  LIBRARIES_TO_LINK
    ${libntn-cho}
    ${libntn-traffic}
    ${libntn-constellation}
    ${libcore}
    ${libnetwork}
    ${libmobility}
    ${libinternet}
    ${libapplications}
    ${libmmwave}
    ${liblte}
)

build_lib_example(
  NAME ntn-cho-full-constellation
  SOURCE_FILES ntn-cho-full-constellation.cc
  LIBRARIES_TO_LINK
    ${libntn-cho}
    ${libntn-traffic}
    ${libntn-constellation}
    ${libntn-observability}
    ${libcore}
    ${libnetwork}
    ${libmobility}
    ${libinternet}
    ${libapplications}
    ${libmmwave}
    ${liblte}
)

build_lib_example(
  NAME ntn-realistic-mobility-demo
  SOURCE_FILES ntn-realistic-mobility-demo.cc
  LIBRARIES_TO_LINK
    ${libntn-cho}
    ${libntn-traffic}
    ${libntn-constellation}
    ${libcore}
    ${libnetwork}
    ${libmobility}
    ${libinternet}
    ${libapplications}
    ${libmmwave}
    ${liblte}
)

# Real UDP traffic to a ground UE on a real mmwave NR cell, handed over by the
# TTE-aware CHO on MEASURED serving SINR.
build_lib_example(
  NAME ntn-cho-handover-traffic
  SOURCE_FILES ntn-cho-handover-traffic.cc
  LIBRARIES_TO_LINK
    ${libntn-cho}
    ${libntn-traffic}
    ${libntn-constellation}
    ${libcore}
    ${libnetwork}
    ${libmobility}
    ${libinternet}
    ${libapplications}
    ${libmmwave}
    ${liblte}
)

# Real-stack flagship: TTE-aware CHO decides on MEASURED mmwave SINR.
build_lib_example(
  NAME ntn-cho-real-stack
  SOURCE_FILES ntn-cho-real-stack.cc
  LIBRARIES_TO_LINK
    ${libntn-cho}
    ${libntn-traffic}
    ${libntn-constellation}
    ${libntn-observability}
    ${libcore}
    ${libnetwork}
    ${libmobility}
    ${libinternet}
    ${libapplications}
    ${libmmwave}
    ${liblte}
)

# The campaign's real-shell mode (--tleTruthFile/--tleNetFile) anchors every
# satellite loaded from real element sets to one simulation epoch through
# Sgp4MobilityModel::SetSimulationEpochUnix, which ntn-constellation gained in
# ns3-ntn-toolkit 8d395304a. On an older ntn-constellation apply
# patches/ntn-constellation-sgp4-simulation-epoch.patch first (see INSTALL.md).
set(ntn_cho_sgp4_header
    ${CMAKE_SOURCE_DIR}/contrib/ntn-constellation/model/sgp4-mobility-model.h)
set(ntn_cho_has_sim_epoch FALSE)
if(EXISTS ${ntn_cho_sgp4_header})
  file(STRINGS ${ntn_cho_sgp4_header} ntn_cho_sim_epoch_hits
       REGEX "SetSimulationEpochUnix")
  if(ntn_cho_sim_epoch_hits)
    set(ntn_cho_has_sim_epoch TRUE)
  endif()
endif()
if(ntn_cho_has_sim_epoch)
  build_lib_example(
    NAME ntn-cho-tte-campaign
    SOURCE_FILES ntn-cho-tte-campaign.cc
    LIBRARIES_TO_LINK
      ${libntn-cho}
      ${libntn-constellation}
      ${libpropagation}
      ${libcore}
      ${libnetwork}
      ${libmobility}
  )
else()
  message(WARNING
    "ntn-cho: contrib/ntn-constellation lacks Sgp4MobilityModel::SetSimulationEpochUnix; "
    "ntn-cho-tte-campaign is not built. Apply "
    "contrib/ntn-cho/patches/ntn-constellation-sgp4-simulation-epoch.patch from the "
    "ns-3 root (git apply) and reconfigure.")
endif()
