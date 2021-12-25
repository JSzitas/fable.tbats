# test script

remove(list=ls())

pkgload::load_all()

lynx_model <- fabletools::model( tsibbledata::pelt,
                                 TBATS(Lynx),
                                 fable::ETS(Lynx)
                                 )

fcst <- forecast( lynx_model, h = 10 )

# cmp <- components(lynx_model)

