# ==============================================================================
# Copyright (c) 2026 IsoMatrix Suite. All rights reserved.
# 
# IsoDecipher is dual-licensed:
# 1. For academic and non-commercial use, it is licensed under the AGPLv3.
# 2. For commercial and enterprise use, a Commercial License is required.
# 
# See the LICENSE file in the project root for more details.
# ==============================================================================


import pandas as pd

# Load the master panel
panel = pd.read_csv("results/panel_features.csv")

# Generate the unique feature names used by IsoDecipher
var_info = panel.copy()
var_info['feature'] = var_info.apply(lambda r: f"{r.gene}_G{r.polyA_group}_{r.user_label}", axis=1)

# Save as the master metadata file
var_info.set_index('feature').to_csv("results/feature_metadata.csv")
print("✅ Master feature metadata saved to results/feature_metadata.csv")