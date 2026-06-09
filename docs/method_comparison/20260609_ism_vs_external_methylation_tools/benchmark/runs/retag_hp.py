import pysam, sys
inbam, outbam, mode = sys.argv[1], sys.argv[2], sys.argv[3]
bf = pysam.AlignmentFile(inbam, 'rb'); of = pysam.AlignmentFile(outbam, 'wb', template=bf)
kept = dropped = 0
for r in bf:
    try: hp = str(r.get_tag('HP'))
    except KeyError: hp = None
    new = None
    if mode == 'somatic':     # ISM HP-axis: HP1(germline) vs HP1-1(somatic)
        if hp == '1': new = 1
        elif hp == '1-1': new = 2
    elif mode == 'family':    # germline-parental: HP1-family vs HP2-family
        if hp in ('1', '1-1'): new = 1
        elif hp in ('2', '2-1'): new = 2
    if new is None:
        dropped += 1; continue
    r.set_tag('HP', new, value_type='i'); of.write(r); kept += 1
bf.close(); of.close()
print(f'mode={mode} kept={kept} dropped={dropped}')
