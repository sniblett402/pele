import argparse
from pele.storage.database import Database
from pele.utils.optim_compatibility import WritePathsampleDB

def main():
    parser = argparse.ArgumentParser(description="Convert a pele database for an angle-axis system into PATHSAMPLE format. The following types of system are currently implemented: 'tetrahedron','otp'")

    parser.add_argument("database", type=str, help="Database file name")
    parser.add_argument("system", type=str, help="System name")
    parser.add_argument("nmol", type=int, help="Number of rigid molecules in the system")

    args = parser.parse_args()

    db = Database(db=args.database, createdb=False)

    if args.system=="tetrahedron":
        import playground
        system = playground.plate_folding.geometric_folding.PlateFolder(args.nmol)
    elif args.system=="otp":
        from pele.pele.angleaxis._otp_cluster import OTPCluster
        system = OTPCluster(args.nmol)
    elif args.system=="tip4p":
        import playground
        system = playground.gmin_tip4p.tip4p_system.TIP4pSystem()
    else:
        raise AttributeError("Unknown system type specified.")

    topology = system.aatopology
    writer = WritePathsampleDB(db,aatopology=topology)
    writer.write_db()


if __name__ == "__main__":
    main()
