import os
from Moltherm.break_sdf import ProductSplit

test_dir = os.path.dirname(os.path.abspath(__file__))

ps = ProductSplit(input_path=test_dir,
                  input_file="cyclohexene.sdf",
                  output_path=test_dir,
                  breaking_bond_sites=[6, 16])
ps.split()
