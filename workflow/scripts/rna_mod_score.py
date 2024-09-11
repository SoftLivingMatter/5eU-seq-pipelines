import pandas as pd
from pathlib import Path
import argparse


class RNAModDataset():
    def __init__(
            self,
            sample_files,
            position_name="Position",
            limit_to_sites=None,
            remove_from_score=None,
        ):
        self.position_name = position_name
        self.remove_from_score = (
            set(int(pos) for pos in remove_from_score)
            if remove_from_score else None
        )
        # all T positions in rDNA
        self.limit_to_sites = (
            set(int(pos) for pos in limit_to_sites)
            if limit_to_sites else None
        )

        self.data = {
            Path(file.name).name: self._read_bedgraph(file)
            for file in sample_files
        }

    def _read_bedgraph(self, file):
        data = pd.read_csv(file,
                           delimiter='\t',
                           names=['samp', 'position', 'end', 'counts'],
                           usecols=[1, 3],
                           )
        if self.limit_to_sites:
            return data[data['position'].isin(self.limit_to_sites)].reset_index(drop=True)
        return data

    def query_positions(self, positions):
        weights = [0.9, 1, 0, 1, 0.9]
        empty = {
                    'fileName': None,
                    self.position_name: None,
                    'ScoreB': None,
                    'ScoreC': None,
                }

        def query_file(file, data, position):
            row = data.index[data['position'] == position]
            if len(row) != 1:
                return empty
            row = row[0]

            # when looking at a few sites, want +/- 2 indices
            if self.limit_to_sites:
                # passing negative values to iloc returns empty df
                start = max(0, row-2)

                # when row is too close to the start of the df, clip start indices
                start_weight = 0
                if row < 2:
                    start_weight += 2 - row

                # when row is too close to the end of the df, clip end indices
                end_weight = 5
                if row > len(data) - 3:
                    end_weight -= 3 - (len(data) - row)

                sub_data = data.iloc[start:row+3].assign(weight=weights[start_weight:end_weight])

            # with the entire dataset, want to look at +/- 2 positions
            else:
                sub_data = data[
                    (position - 2 <= data['position']) &
                    (data['position'] <= position + 2)
                ].copy()

                if position not in sub_data['position'].values:
                    return empty

                # map weight to proper index based on position
                sub_data['weight'] = [weights[2+abs(position-p)]
                    for p in sub_data['position']]
            
            # remove_from_score sites are removed from calculation
            if self.remove_from_score:
                sub_data = sub_data[~sub_data['position'].isin(self.remove_from_score)]

            if sub_data['weight'].sum() == 0:
                return empty

            weighted_average = (sub_data['counts'] * sub_data['weight']).sum() / sub_data['weight'].sum()
            return {
                'fileName': file,
                self.position_name: position,
                'ScoreB': abs(data.iloc[row]['counts'] - weighted_average) / (data.iloc[row]['counts'] + 1),
                'ScoreC': 1 - (data.iloc[row]['counts'] / (weighted_average)),
            }
            
        return pd.DataFrame([query_file(file, data, position)
            for position in positions
            for file, data in self.data.items()
        ]).dropna().astype({self.position_name: int})  # remove positions that aren't found

def parse():
    parser = argparse.ArgumentParser(
        prog='mod_score',
        description='Caculate scoreB and scoreC for methyl experiments',
    )
    parser.add_argument('--query-sites', type=argparse.FileType('r'), required=True)
    parser.add_argument('--limit-sites', type=argparse.FileType('r'))
    parser.add_argument('--exclude-score', type=argparse.FileType('r'))
    parser.add_argument('--column-name', default='Position')
    parser.add_argument('--output', required=True, type=argparse.FileType('w'))
    parser.add_argument('sample_files', type=argparse.FileType('r'), nargs='+')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    dataset = RNAModDataset(
        args.sample_files,
        position_name=args.column_name,
        remove_from_score=args.exclude_score,
        limit_to_sites=args.limit_sites,
    )
    result = dataset.query_positions(int(line) for line in args.query_sites)
    result.to_csv(args.output, index=False)
