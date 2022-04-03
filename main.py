import sys
from src.pipeline import Pipeline

if __name__ == "__main__":
    print("Decouphage 0.0.1")
    input_path = sys.argv[1]
    Pipeline(input_path)
