import numpy as np
import shapely

from threedigrid_builder.base import Array, search

__all__ = ["ExchangeLines"]


class ExchangeLine:
    id: int
    the_geom: shapely.Geometry
    channel_id: int
    exchange_level: float


class ExchangeLines(Array[ExchangeLine]):
    @property
    def channel_mapping(self) -> "ExchangeLineChannels":
        """Return 2 exchange lines ids for each channel id"""
        if not hasattr(self, "_channel_mapping"):
            channel_ids = np.unique(self.channel_id)
            result = np.full((len(channel_ids), 2), fill_value=-9999, dtype=np.int32)
            for i, mask in enumerate(self.split_in_two(self.channel_id)):
                result[:, i] = self.index_to_id(
                    search(
                        self.channel_id,
                        channel_ids,
                        mask=mask,
                        assume_ordered=False,
                        check_exists=False,
                    )
                )
            self._channel_mapping = ExchangeLineChannels(
                id=channel_ids,
                exchange_line_id=result[:, 0],
                secondary_exchange_line_id=result[:, 1],
            )
        return self._channel_mapping

    def get_for_channel_id(self, channel_id, is_primary):
        if len(self) == 0:
            return np.full_like(channel_id, fill_value=-9999)
        index = self.channel_mapping.id_to_index(channel_id)
        if is_primary:
            ids = self.channel_mapping.exchange_line_id
        else:
            ids = self.channel_mapping.secondary_exchange_line_id
        return np.where(index != -9999, np.take(ids, index, mode="clip"), -9999)


class ExchangeLineChannel:
    id: int
    exchange_line_id: int
    secondary_exchange_line_id: int


class ExchangeLineChannels(Array[ExchangeLineChannel]):
    pass
