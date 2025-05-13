energy_loss_data = [
        ("5%", "bifacial", -0.560, 0.133, 40.2, 0.70, -0.00268, 0.361),
        ("5%", "monofacial", -0.550, 0.138, 43.4, 0.71, -0.00282, 0.388),
        ("10%", "bifacial",  -0.485, 0.171, 46.2, 0.72, -0.00437, 0.575),
        ("10%", "monofacial", -0.441, 0.198, 48.7, 0.74, -0.00476, 0.621),
        ("15%", "bifacial",  -0.414, 0.207, 49.9, 0.74, -0.00576, 0.762),
        ("15%", "monofacial", -0.371, 0.208, 51.5, 0.75, -0.00633, 0.825)
    ]
    #creating a dataframe from the data in energy_loss_data
    energy_loss_df = DataFrame(
        energy_loss = [row[1] for row in energy_loss_data],
        solar_type = [row[2] for row in energy_loss_data],
        P = [row[3] for row in energy_loss_data],
        k = [row[4] for row in energy_loss_data],
        α_0 = [row[5] for row in energy_loss_data],
        GCR_0 = [row[6] for row in energy_loss_data],
        m = [row[7] for row in energy_loss_data],
        b = [row[8] for row in energy_loss_data] 
    )
    println(energy_loss_df)