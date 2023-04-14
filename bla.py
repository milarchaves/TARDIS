count = 0
with open ('ogee.txt') as ogee:
    ogee = ogee.readlines()
    with open ('tardis.txt') as tardis:
        for line in tardis:
            if line in ogee:
                count += 1
                print(line)
print(count)