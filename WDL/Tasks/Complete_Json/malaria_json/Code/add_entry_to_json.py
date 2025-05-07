import json
import sys

def add_column_to_json(txt_file, json_file):
    try:
        with open(txt_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: File not found: {txt_file}")
        sys.exit(1)

    if len(lines) < 2:
        print("Error: TXT file must contain a header and at least one value.")
        sys.exit(1)

    key = lines[0]
    values = lines[1:]

    value_to_store = values[0] if len(values) == 1 else values

    # Load existing JSON or create new
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        data = {}

    data[key] = value_to_store

    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)

    print(f"Added '{key}': {value_to_store} to {json_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_column_to_json.py <txt_file> <json_file>")
        sys.exit(1)

    txt_file = sys.argv[1]
    json_file = sys.argv[2]
    add_column_to_json(txt_file, json_file)

