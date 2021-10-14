#!/usr/bin/env python3

# file://check_mast.py

# Kenneth Mighell
# SETI Institute

# =============================================================================


def check_mkpy3():
    """
Purpose: check MAST availibilty with mkpy3_bad_radec_bug_v1.py

Kenneth Mighell
SETI Institute
    """
    import subprocess
    import datetime
    import time

    # ASCII color codes
    CEND = "\33[0m"
    # CRED = '\33[31m'  # red
    CRED = "\33[1;31m"  # bold red
    # CGREEN = '\33[32m'  # green
    CGREEN = "\33[1;32m"  # bold green

    file = "./mkpy3_bad_radec_bug_v1.py"
    print(file)

    good = 0
    bad = 0
    nloops = 10
    for k in range(nloops):
        timestamp = str(datetime.datetime.now())
        print(k, timestamp, "sent...")
        result = subprocess.run(["python3", file], capture_output=True)
        timestamp = str(datetime.datetime.now())
        print(k, timestamp, "...received")
        returncode = result.returncode
        if returncode != 0:
            bad += 1
            print(CRED + "***** ERROR FOUND *****" + CEND)
            print()
        else:
            good += 1
            print(CGREEN + "OK" + CEND)
            # fi
        timestamp = str(datetime.datetime.now())
        # print(k, file, timestamp)
        delay_sec = 10
        if k < (nloops - 1):
            time.sleep(delay_sec)
            print("wait", delay_sec, "seconds...")
        # rof

    print()
    print(good, "=good")
    print(bad, "=bad")


# =============================================================================


if __name__ == "__main__":
    check_mkpy3()
    # fi

# EOF
