list(APPEND multiAttrSort-core_sources
    composer.cpp
    chain_composer.cpp
    stitch_composer.cpp	
    plan_enum_rrs.cpp
    )

add_library(multiAttrSort-core STATIC ${multiAttrSort-core_sources})
