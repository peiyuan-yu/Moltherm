import os
from Moltherm.break_sdf import ProductSplit

test_dir = os.path.dirname(os.path.abspath(__file__))

ps = ProductSplit.from_file(input_path=test_dir,
                            input_file="cyclohexene.sdf",
                            output_path=test_dir)
ps.split(break_bond_sites=[5, 15])
