import ase
from ase.calculators.orca import ORCA, OrcaProfile
import argparse
from ase.io import read, write
import os
from ase.optimize import LBFGS

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="path to input configs", required=True)
    parser.add_argument("--output", help="path to save output configs to", required=True)
    parser.add_argument("--charge", help="charge", default=0, type=int)
    parser.add_argument("--mult", help="spin multiplicity", default=1, type=int)
    parser.add_argument("--task", help="task", required=True)
    parser.add_argument('--id', help='id', default=0)
    return parser.parse_args()

def get_calc(orca_path, charge=0, mult=1, job_id=0):
    calc = ORCA(
        profile=OrcaProfile(orca_path),
        label=os.path.join("orca"),
        charge=charge, mult=mult,task='gradient',
        orcasimpleinput='PBE D3 def2-TZVP TightSCF',
        orcablocks='%pal nprocs 10 end \n %scf ConvForced=1 end'
    )
    return calc

def single_point(orca_path, config, charge, mult, job_id, prefix="DFT_"):
    calc = get_calc(orca_path, charge, mult, job_id)
    calc.calculate(atoms=config, properties={"energy"}, system_changes=None)
    print(calc.results)
    #os.system("rm orca_property.txt orca.*")
    config.info.update(
        {
            f"{prefix}energy": calc.results["energy"]
        }
    )
    return config


def optimize(config, charge, mult, out, job_id,fmax=0.01, nsteps=100):
    calc = get_calc(charge, mult)
    config.calc = calc
    opt = LBFGS(config)
    opt.run(fmax=fmax, steps=nsteps)
    write(out, opt.atoms)

def main():
    args = parse_args()
    config = read(args.input, args.id)
    if args.task == "single_point":
        config = single_point('/rds/user/pcpt3/hpc-work/orca/orca', config, args.charge, args.mult, args.id)
        write(args.output, config, append=True)
        #os.system("rm orca.* orca_property.txt")
    elif args.task == "optimize":
        optimize(
            configs[0], 
            charge=args.charge, 
            mult=args.mult,
            out=args.output,
            job_id = args.id
        )

if __name__ == "__main__":
    main()

