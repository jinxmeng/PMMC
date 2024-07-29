#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-05-24, 14:45:51
# modified date: 2024-05-24, 14:45:51

import sys, sqlite3
import argparse

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", metavar="in_f", type=str, help="in_f sqlite3-format db")
    parser.add_argument("--tl", action="store_true", help="out all tables list in db")
    parser.add_argument("-tn", metavar="table_name", type=str, default="none", help="out a table with specify names")
    parser.add_argument("-o", metavar="out_f", type=str, help="out a table names")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args

def get_tl(i):
    conn = sqlite3.connect(i) # 连接到 SQLite 数据库
    cursor = conn.cursor()  # 创建一个游标对象
    # 查看数据库中的所有表名
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    print("Tables in the database:")
    for table in tables:
        print(table[0])
    # 关闭游标和连接
    cursor.close()
    conn.close()


def out_table(i, tn, o):
    conn = sqlite3.connect(i) # 连接到 SQLite 数据库
    cursor = conn.cursor()  # 创建一个游标对象
    out_f = open(o, "w")
    # 打印表的列名（可选）
    cursor.execute(f"PRAGMA table_info({tn});")
    columns = [column[1] for column in cursor.fetchall()]
    # print(f"\nColumns in table '{tn}': {columns}")
    out_f.write("\t".join(columns) + "\n")
    # 获取所有行            
    cursor.execute(f"SELECT * FROM {tn}")
    rows = cursor.fetchall()
    for row in rows:
        out_f.write("\t".join([str(x) for x in row]) + "\n")
    # 关闭游标和连接
    cursor.close()
    conn.close()


def main(i, tl, tn, o):
    if tl is True:
        get_tl(i)
        sys.exit(0)

    if tn == "none":
        sys.exit("Do nothing!")
    else:
        out_table(i, tn, o)

if __name__ == "__main__":
    args = get_args()
    main(args.i, args.tl, args.tn, args.o)

