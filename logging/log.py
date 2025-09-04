import serial
import csv

# Adjust device and baudrate to your setup
ser = serial.Serial("/dev/ttyUSB0", 115200, timeout=1)

with open("timings.csv", "w", newline="") as f:
    writer = csv.writer(f)

    try:
        print("Logging... Press Ctrl+C to stop.")
        while True:
            line = ser.readline().decode(errors="ignore").strip()
            if not line:
                continue

            # Expect lines like: 3125.432,50.123,1840.778,12123.111
            parts = line.split(",")
            if len(parts) == 4:   # expected 4 columns
                writer.writerow(parts)
                print(parts)  # echo to console
            else:
                print(line)
                continue
    except KeyboardInterrupt:
        print("\nStopped logging.")
