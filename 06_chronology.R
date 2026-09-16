library(serac)

#serac adm TL09_22_01----

TL09_22_01_adm <- serac(name="TL09_22_01",
                      model = c("CFCS", "CRS"),
                      coring_yr = 2022,
                      plot_Pb=TRUE,
                      plot_Cs=TRUE,
                      NWT = c(120, 180),
                      plotpdf = TRUE,
                      inst_deposit = c(90, 100, 225, 345),
                      Hemisphere=c("NH"),
                      mass_depth = FALSE,
                      plotphoto = TRUE,
                      minphoto = 0,
                      maxphoto = 835,
                      plot_CFCS_regression = FALSE,
                      plot_unit = "cm",
                      mycex = 1.4
                      )


#serac_input_formatting("TL09_B_60")


tl09_B_60_adm <- serac("TL09_B_60",
                        model = c("CFCS", "CRS"),
                        coring_yr = 2022,
                       plotphoto = TRUE,
                       minphoto = 0,
                       maxphoto = 770,
                       plot_Pb=TRUE,
                        plot_Cs=TRUE,
                        NWT = c(60, 120),
                        plotpdf = TRUE,
                       plot_CFCS_regression = FALSE,
                        inst_deposit = c(160, 200, 300, 460),
                        Hemisphere=c("NH"),
                        mass_depth = FALSE,
                       plot_unit = "cm",
                       mycex = 1.4
)
