from artemis.io.argparser import parse
import sys

def main():

    args = parse()

    args.func(args)


if __name__ == "__main__":
    main()
    sys.exit(0)
