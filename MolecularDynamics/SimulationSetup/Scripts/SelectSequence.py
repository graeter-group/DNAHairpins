

def print_sequence(top_seq_in, bot_seq_in, current_sample_in):
    unpaired = current_sample_in[0]
    stem = current_sample_in[1]
    loop = current_sample_in[2]

    mid_point_top = int(len(top_seq_in)/2)
    left_ds = top_seq_in[mid_point_top-(stem+int(loop/2))-50:mid_point_top-(stem+int(loop/2))]
    left_hp = top_seq_in[mid_point_top-(stem+int(loop/2)):mid_point_top-int(loop/2)]
    loop_hp = top_seq_in[mid_point_top-int(loop/2):mid_point_top+int(loop/2)]
    right_hp = top_seq_in[mid_point_top+int(loop/2):mid_point_top+stem+int(loop/2)]
    right_ds = top_seq_in[mid_point_top+(stem+int(loop/2)): mid_point_top+(stem+int(loop/2))+50:]
    mid_point_bot = len(bot_seq_in)/2

    right_ds_bot = bot_seq_in[int(mid_point_bot-unpaired/2)-50: int(mid_point_bot-unpaired/2)][::-1]
    left_ds_bot = bot_seq_in[int(mid_point_bot+unpaired/2): int(mid_point_bot+unpaired/2)+50][::-1]

    unpaired_bot = bot_seq_in[int(mid_point_bot-unpaired/2): int(mid_point_bot+unpaired/2)]

    print(f"HP {unpaired} {stem} {loop}")
    print(f"TOP {left_ds} {left_hp} {loop_hp} {right_hp} {right_ds}")
    print(f"BOT {left_ds_bot} {unpaired_bot} {right_ds_bot}")
    print("\n")


raw_sequence_file = open("RawSequences.txt")

current_sample = None
top_seq = ""
bot_seq = ""

for line in raw_sequence_file:
    if line.startswith("HP"):
        line_data = [int(a) for a in line.strip("\n").split(" ")[1:]]
        if current_sample and line_data != current_sample:
            print_sequence(top_seq, bot_seq, current_sample)

        current_sample = line_data
        continue
    if line.startswith("TOP"):
        top_seq = line.strip("\n").split(" ")[1]
        continue
    if line.startswith("BOT"):
        bot_seq = line.strip("\n").split(" ")[1]
        continue

print_sequence(top_seq, bot_seq, current_sample)
