

N_DESIRED_BIRTHS = 3
N_IVF_PER_COUPLE = 2  # Number of IVF cycles per couple
N_EMBRYOS_PER_IVF = 6  # Number of viable embryos per IVF cycle
FRACTION_INFERTILE = 0.05  # Fraction of couples that are infertile
BIRTH_RATE = 0.65  # Live birth rate per transferred viable embryo
RANDOM_STATE = 0  # Random seed for reproducibility


N_GENERATIONS = 8


DEGREE_UNRELATED = 4
# child would be offspring of... (parents would have... )
# 1: siblings (shared parents)
# 2: first cousins (shared grandparents)
# 3: second cousins (shared great grandparents)
# 4: third cousins (shared great-great grandparents) - unrelated



# The below comes from running: grep -P "^chr[1-9XY]" Homo_sapiens_assembly38.fasta.fai | cut -f1,2
CHROMOSOME_LENGTHS = {
    "1": 248956422, 
    "2": 242193529, 
    "3": 198295559, 
    "4": 190214555, 
    "5": 181538259, 
    "6": 170805979, 
    "7": 159345973, 
    "8": 145138636, 
    "9": 138394717, 
    "10": 133797422, 
    "11": 135086622, 
    "12": 133275309, 
    "13": 114364328, 
    "14": 107043718, 
    "15": 101991189, 
    "16": 90338345, 
    "17": 83257441, 
    "18": 80373285, 
    "19": 58617616, 
    "20": 64444167, 
    "21": 46709983, 
    "22": 50818468, 
    "X": 156040895, 
    "Y": 57227415
}






