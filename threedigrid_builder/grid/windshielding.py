from threedigrid_builder.base import array_of
from threedigrid_builder.base import search
from threedigrid_builder.constants import ContentType

import numpy as np


__all__ = ["Windshieldings"]


class Windshielding:
    id: int
    channel_id: int
    north: float
    northeast: float
    east: float
    southeast: float
    south: float
    southwest: float
    west: float
    northwest: float


@array_of(Windshielding)
class Windshieldings:
    def apply_to_lines(self, lines, channels):
        """Set windshielding factors on lines. Only channels can have windshielding,
        as only channels can have wind applied to them

            Args:
                lines (Lines)
                channels (Channels)
        """

        if len(channels.id) == 0 or len(self.id) == 0:
            return

        if not np.in1d(channels.id, self.channel_id).all():
            missing = set(self.channel_id) - set(channels.id)
            raise ValueError(f"Windshielding given for non existing channels {missing}")

        mask = lines.content_type == ContentType.TYPE_V2_CHANNEL
        sorter = np.argsort(self.channel_id)
        chn_idx = np.digitize(
            lines.content_pk[mask], self.channel_id[sorter], right=True
        )

        lines.windshieldings[:, 0][mask] = self.north[chn_idx]
        lines.windshieldings[:, 1][mask] = self.northeast[chn_idx]
        lines.windshieldings[:, 2][mask] = self.east[chn_idx]
        lines.windshieldings[:, 3][mask] = self.southeast[chn_idx]
        lines.windshieldings[:, 4][mask] = self.south[chn_idx]
        lines.windshieldings[:, 5][mask] = self.southwest[chn_idx]
        lines.windshieldings[:, 6][mask] = self.west[chn_idx]
        lines.windshieldings[:, 7][mask] = self.northwest[chn_idx]
